"""Exercise draft preservation and save races without browser dependencies."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CLIENT = ROOT / "scripts/publication/candidate_review_notes.js"
NODE = shutil.which("node")

HARNESS = r"""
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const source = fs.readFileSync(process.argv[1], "utf8");
const tick = async () => { for (let i = 0; i < 5; i++) await new Promise(setImmediate); };
const deferred = () => {
  let resolve;
  const promise = new Promise(done => { resolve = done; });
  return { promise, resolve };
};
const response = (body, status = 200) => ({
  ok: status >= 200 && status < 300, status, json: async () => body,
});
const entry = (text, revision) => ({text, revision, title: "Figure 5", updated_at: "now"});
const notes = value => response({schema_version: 1, notes: value ? {"figure-5": value} : {}});
const key = "desmognathus.figure-review.notes.v1:/:figure-5";

function page(handler, storage = new Map(), storageBlocked = false) {
  const timers = new Map();
  let timerId = 0;
  class Element {
    constructor() {
      this.listeners = {}; this.value = ""; this.textContent = ""; this.disabled = true;
      this.hidden = true;
    }
    addEventListener(name, fn) { (this.listeners[name] ||= []).push(fn); }
    emit(name, event = {}) { return Promise.all((this.listeners[name] || []).map(fn => fn(event))); }
    focus() { this.focused = true; }
    closest() { return section; }
    set innerHTML(_) { throw new Error("Do not insert note content as HTML"); }
    click() { this.clicked = true; return this.emit("click"); }
    remove() { this.removed = true; }
  }
  const input = new Element(); input.dataset = {noteId: "figure-5"};
  const heading = new Element(); heading.textContent = "Figure 5 · Example";
  const status = new Element();
  const conflict = new Element();
  const remote = new Element();
  const focus = new Element();
  const useSaved = new Element();
  const keepDraft = new Element();
  const summary = new Element();
  const save = new Element();
  const exported = new Element();
  const section = {querySelector: selector => ({
    "h2": heading, ".note-status": status, ".note-conflict": conflict, ".remote-note": remote,
    '[data-note-action="focus"]': focus,
    '[data-note-action="use-saved"]': useSaved,
    '[data-note-action="keep-draft"]': keepDraft,
  })[selector]};
  const calls = [];
  const window = new Element();
  const downloads = [];
  let blob;
  class BrowserURL extends URL {
    static createObjectURL(value) { blob = value; return "blob:test"; }
    static revokeObjectURL() {}
  }
  const context = {
    document: {
      readyState: "complete",
      querySelector: selector => ({
        "#notes-summary": summary, "#save-notes": save, "#export-notes": exported,
      })[selector],
      querySelectorAll: () => [input],
      createElement: () => new Element(),
      body: {appendChild: value => downloads.push(value)},
    },
    window, location: {pathname: "/", href: "http://127.0.0.1:8766/"},
    localStorage: {
      getItem: name => { if (storageBlocked) throw new Error("blocked"); return storage.get(name) || null; },
      setItem: (name, value) => { if (storageBlocked) throw new Error("blocked"); storage.set(name, value); },
      removeItem: name => { if (storageBlocked) throw new Error("blocked"); storage.delete(name); },
    },
    fetch: (url, options) => { calls.push({url, options}); return handler(url, options); },
    setTimeout: (fn, ms) => { const id = ++timerId; timers.set(id, {fn, ms}); return id; },
    clearTimeout: id => timers.delete(id),
    AbortController, URL: BrowserURL, Blob,
  };
  vm.runInNewContext(source, context);
  return {
    input, status, conflict, remote, useSaved, keepDraft, summary, save, exported,
    calls, storage, window, downloads, blob: () => blob,
    type: text => { input.value = text; return input.emit("input"); },
    runTimers: ms => {
      for (const [id, timer] of [...timers]) {
        if (timer.ms === ms) { timers.delete(id); timer.fn(); }
      }
    },
  };
}
"""


def run_client(scenario: str) -> None:
    if NODE is None:
        pytest.skip("Node is required for the browser notes state checks")
    result = subprocess.run(
        [NODE, "-e", HARNESS + "\n(async () => {\n" + scenario + "\n})().catch(error => {\n"
         "console.error(error); process.exitCode = 1;\n});", str(CLIENT)],
        check=False,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_edits_during_request_serialize_with_the_acknowledged_revision() -> None:
    run_client(r"""
      const first = deferred(); const second = deferred();
      let puts = 0;
      const app = page((url, options) => options.method === "PUT"
        ? (++puts === 1 ? first.promise : second.promise) : notes());
      await tick();
      assert.equal(app.input.disabled, false);
      await app.type("First draft");
      app.input.emit("blur");
      await tick();
      await app.type("Second draft typed before acknowledgement");
      assert.equal(puts, 1);
      first.resolve(response({entry: entry("First draft", 1)}));
      await tick();
      assert.equal(puts, 2);
      const bodies = app.calls.filter(call => call.options.method === "PUT")
        .map(call => JSON.parse(call.options.body));
      assert.deepEqual(bodies.map(body => body.base_revision), [0, 1]);
      assert.equal(bodies[1].text, "Second draft typed before acknowledgement");
      assert.match(app.status.textContent, /Saving to project/);
      assert.equal(JSON.parse(app.storage.get(key)).text, bodies[1].text);
      second.resolve(response({entry: entry(bodies[1].text, 2)}));
      await tick();
      assert.equal(app.status.textContent, "Saved to project");
      assert.equal(app.storage.has(key), false);
    """)


def test_offline_browser_draft_survives_refresh_and_explicit_retry() -> None:
    run_client(r"""
      const storage = new Map();
      const offline = page(() => Promise.reject(new Error("offline")), storage);
      await tick();
      await offline.type("Check axis labels <script>alert(1)</script>");
      assert.match(offline.status.textContent, /Draft saved in this browser/);
      let prevented = false;
      await offline.window.emit("beforeunload", {preventDefault: () => { prevented = true; }});
      assert.equal(prevented, true);
      let online = false;
      const refreshed = page((url, options) => {
        if (!online) return Promise.reject(new Error("offline"));
        if (options.method !== "PUT") return notes();
        return response({entry: entry(JSON.parse(options.body).text, 1)});
      }, storage);
      await tick();
      assert.equal(refreshed.input.value, offline.input.value);
      online = true;
      await refreshed.save.click();
      await tick();
      assert.equal(refreshed.status.textContent, "Saved to project");
      assert.equal(storage.has(key), false);
      assert.equal(refreshed.calls.filter(call => call.options.method === "PUT").length, 1);
    """)


def test_revision_conflict_preserves_draft_and_requires_an_explicit_choice() -> None:
    run_client(r"""
      const storage = new Map([[key, JSON.stringify({text: "My local draft", base_revision: 1})]]);
      const app = page((url, options) => options.method === "PUT"
        ? response({entry: entry(JSON.parse(options.body).text, 4)})
        : notes(entry("Other saved draft", 3)), storage);
      await tick();
      assert.equal(app.conflict.hidden, false);
      assert.equal(app.input.value, "My local draft");
      assert.equal(app.remote.textContent, "Other saved draft");
      await app.save.click();
      assert.equal(app.calls.filter(call => call.options.method === "PUT").length, 0);
      await app.keepDraft.click();
      await tick();
      assert.equal(app.conflict.hidden, true);
      const request = app.calls.find(call => call.options.method === "PUT");
      assert.equal(JSON.parse(request.options.body).base_revision, 3);
      assert.equal(app.status.textContent, "Saved to project");
    """)


def test_lost_acknowledgement_does_not_mark_a_reverted_draft_as_saved() -> None:
    run_client(r"""
      let puts = 0;
      const app = page((url, options) => {
        if (options.method !== "PUT") return notes(entry("Original", 1));
        if (++puts === 1) return Promise.reject(new Error("Acknowledgement lost"));
        return response({error: "conflict", entry: entry("Changed on server", 2)}, 409);
      });
      await tick();
      await app.type("Changed on server");
      await app.input.emit("blur");
      await app.type("Original");
      assert.match(app.status.textContent, /Draft saved in this browser/);
      assert.equal(app.storage.has(key), true);
      await app.save.click();
      assert.equal(app.conflict.hidden, false);
      assert.equal(app.input.value, "Original");
      assert.equal(app.remote.textContent, "Changed on server");
    """)


def test_conflict_response_and_use_saved_preserve_an_empty_project_note() -> None:
    run_client(r"""
      const app = page((url, options) => options.method === "PUT"
        ? response({error: "conflict", entry: entry("", 5)}, 409)
        : notes(entry("Original saved note", 3)));
      await tick();
      await app.type("My changed draft");
      await app.input.emit("blur");
      assert.equal(app.conflict.hidden, false);
      assert.equal(app.input.value, "My changed draft");
      assert.match(app.remote.textContent, /project note is empty/);
      await app.useSaved.click();
      assert.equal(app.input.value, "");
      assert.equal(app.conflict.hidden, true);
      assert.equal(app.storage.has(key), false);
      assert.equal(app.status.textContent, "Saved to project");
    """)


def test_removed_project_note_can_be_resolved_and_recreated() -> None:
    run_client(r"""
      let puts = 0;
      const app = page((url, options) => {
        if (options.method !== "PUT") return notes(entry("Previous saved note", 3));
        if (++puts === 1) return response({error: "conflict", entry: null}, 409);
        return response({entry: entry(JSON.parse(options.body).text, 1)});
      });
      await tick();
      await app.type("My changed draft");
      await app.input.emit("blur");
      assert.equal(app.conflict.hidden, false);
      assert.equal(app.input.value, "My changed draft");
      await app.keepDraft.click();
      await tick();
      const request = app.calls.filter(call => call.options.method === "PUT")[1];
      assert.equal(JSON.parse(request.options.body).base_revision, 0);
      assert.equal(app.status.textContent, "Saved to project");
    """)


def test_unavailable_browser_storage_is_reported_and_export_keeps_the_draft() -> None:
    run_client(r"""
      const app = page(() => Promise.reject(new Error("offline")), new Map(), true);
      await tick();
      await app.type("Keep this unsaved note\n<script>alert(1)</script>");
      assert.match(app.status.textContent, /Browser backup unavailable/);
      assert.equal(app.exported.disabled, false);
      await app.exported.click();
      const content = await app.blob().text();
      assert.match(content, /## Figure 5 · Example/);
      assert.match(content, /http:\/\/127\.0\.0\.1:8766\/#figure-5/);
      assert.match(content, /Keep this unsaved note\n<script>alert\(1\)<\/script>/);
      assert.match(content, /has not been saved to the project/);
      assert.equal(app.downloads[0].clicked, true);
      assert.equal(app.downloads[0].removed, true);
    """)


def test_saved_notes_on_removed_panels_remain_counted_and_exportable() -> None:
    run_client(r"""
      const removed = {...entry("Remove this panel\nKeep the remaining four.", 2),
        title: "Figure 10b · Development mode"};
      const app = page(() => response({schema_version: 1, notes: {
        "figure-10-b": removed,
        "figure-10-x": {...entry("   ", 1), title: "Empty removed note"},
      }}));
      await tick();
      assert.equal(app.input.value, "");
      assert.equal(app.exported.disabled, false);
      assert.match(app.summary.textContent, /1 figure\/panel has notes/);
      assert.match(app.summary.textContent, /1 note on a removed panel/);
      assert.match(app.summary.textContent, /all notes saved to project/);
      await app.exported.click();
      const content = await app.blob().text();
      assert.match(content, /## Figure 10b · Development mode/);
      assert.match(content, /Saved note\. Panel removed from current gallery\./);
      assert.ok(content.includes(removed.text));
      assert.doesNotMatch(content, /View figure or panel|#figure-10-b|Empty removed note/);
      assert.equal(app.calls.filter(call => call.options.method === "PUT").length, 0);
      assert.equal(app.storage.size, 0);
      let prevented = false;
      await app.window.emit("beforeunload", {preventDefault: () => { prevented = true; }});
      assert.equal(prevented, false);
    """)


def test_removed_notes_export_alongside_current_drafts_without_duplicates() -> None:
    run_client(r"""
      const app = page(() => response({schema_version: 1, notes: {
        "figure-5": entry("Current saved note", 1),
        "figure-10-b": {...entry("Removed development panel note", 2), title: "Development mode"},
        "figure-10-f": {...entry("Another removed panel note", 1), title: "Previous sixth panel"},
      }}));
      await tick();
      await app.type("Updated current draft");
      assert.match(app.summary.textContent, /3 figures\/panels have notes/);
      assert.match(app.summary.textContent, /2 notes on removed panels/);
      assert.match(app.summary.textContent, /1 unsaved/);
      await app.exported.click();
      const content = await app.blob().text();
      assert.equal(content.match(/^## /gm).length, 3);
      assert.match(content, /http:\/\/127\.0\.0\.1:8766\/#figure-5/);
      assert.match(content, /Updated current draft/);
      assert.doesNotMatch(content, /Current saved note|#figure-10-[bf]/);
      assert.equal(content.match(/Panel removed from current gallery/g).length, 2);
      assert.equal(content.match(/has not been saved to the project/g).length, 1);
      assert.equal(JSON.parse(app.storage.get(key)).text, "Updated current draft");
      assert.equal(app.calls.filter(call => call.options.method === "PUT").length, 0);
    """)
