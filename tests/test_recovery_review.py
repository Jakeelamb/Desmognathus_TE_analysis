"""Exercise the shipped review page without a browser or scientific artifacts."""

from __future__ import annotations

import subprocess
from html.parser import HTMLParser
from pathlib import Path

import pytest

PAGE = Path(__file__).resolve().parents[1] / "scripts/microscopy/recovery_review.html"

HARNESS = r"""
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const {webcrypto} = require('node:crypto');
const source = fs.readFileSync(process.argv[1], 'utf8');
const scripts = [...source.matchAll(/<script[^>]*>([\s\S]*?)<\/script>/g)].map(match => match[1]);
const readyChecks = [];
const tick = async () => {
  for (let i = 0; i < 200; i++) {
    await new Promise(resolve => setTimeout(resolve, 1));
    if (readyChecks.every(check => check())) return;
  }
  throw Error('Review page did not finish its asynchronous image/state load');
};
const normalize = value => JSON.parse(JSON.stringify(value));
const row = (id = 'case_01', state = 'a'.repeat(64)) => ({
  case_id: id, raw: `images/${id}/raw.png`, mask: `images/${id}/mask.png`,
  width: 30, height: 20, state_sha256: state,
});
function page(rows = [row(), row('case_02')], storage = new Map(), options = {}) {
  class Element {
    constructor(id = '') {
      this.id = id; this.value = ''; this.checked = false; this.disabled = false;
      this.hidden = false; this.textContent = ''; this.style = {}; this.children = [];
      this.listeners = {}; this.classList = {toggle() {}};
    }
    addEventListener(name, handler) { (this.listeners[name] ||= []).push(handler); }
    emit(name, event = {}) { return Promise.all((this.listeners[name] || []).map(handler => handler(event))); }
    appendChild(child) { this.children.push(child); }
    click() { this.clicked = true; return this.emit('click'); }
    remove() {}
    set innerHTML(_) { throw new Error('Unsafe HTML insertion'); }
  }
  const ids = [...source.matchAll(/\bid="([^"]+)"/g)].map(match => match[1]);
  const elements = Object.fromEntries(ids.map(id => [id, new Element(id)]));
  elements.zoom.value = '1'; elements.opacity.value = '35'; elements['show-mask'].checked = true;
  const decisions = ['keep', 'problem', 'unsure'].map(value => Object.assign(new Element(), {value}));
  const reasons = ['complete_isolated_nucleus', 'incomplete_boundary', 'multiple_nuclei', 'overlap_or_debris', 'edge_truncation', 'uncertain_boundary', 'other'].map(value => Object.assign(new Element(), {value}));
  const document = new Element(); document.body = new Element();
  document.getElementById = id => elements[id];
  document.querySelectorAll = selector => selector.includes('decision') ? decisions : reasons;
  document.createElement = () => new Element();
  const images = [], calls = [], downloads = [];
  class Image {
    set src(value) {
      images.push(value);
      const match = rows.find(item => item.raw === value || item.mask === value);
      this.naturalWidth = match?.width || 30; this.naturalHeight = match?.height || 20;
      if (options.badDimensions) this.naturalWidth++;
      queueMicrotask(() => options.imageFailure ? this.onerror() : this.onload());
    }
  }
  class BrowserURL extends URL {
    static createObjectURL(blob) { downloads.push(blob); return 'blob:review'; }
    static revokeObjectURL() {}
  }
  const context = vm.createContext({
    document, Image, Blob, URL: BrowserURL, crypto: webcrypto, TextEncoder,
    setTimeout: callback => { callback(); },
    fetch: async (url, config) => { calls.push({url, config}); return {ok: !options.fetchFailure, text: async () => JSON.stringify(rows)}; },
    localStorage: {
      getItem: key => { if (options.storageBlocked) throw Error('blocked'); return storage.get(key) || null; },
      setItem: (key, value) => { if (options.storageBlocked) throw Error('blocked'); storage.set(key, value); },
    },
  });
  scripts.forEach(script => vm.runInContext(script, context));
  const run = code => vm.runInContext(code, context);
  readyChecks.push(() => run('imageReady || el("message").textContent.includes("failed") || el("message").textContent.includes("queue could not load")'));
  const choose = async (decision, reason, note = '') => {
    decisions.forEach(input => { input.checked = input.value === decision; });
    reasons.forEach(input => { input.checked = input.value === reason; });
    elements.note.value = note; await elements['review-form'].emit('input');
  };
  return {elements, decisions, reasons, document, images, calls, downloads, storage, run, choose};
}
"""


def run_page(scenario: str) -> None:
    result = subprocess.run(
        [
            "node", "-e",
            HARNESS + "\n(async () => {\n" + scenario
            + "\n})().catch(error => { console.error(error); process.exitCode = 1; });",
            str(PAGE),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_starts_blank_reads_only_blinded_case_list_and_shows_native_pixels():
    run_page(r"""
      const app = page(); await tick();
      assert.equal(app.elements.queue.hidden, false);
      assert.equal(app.decisions.some(input => input.checked), false);
      assert.equal(app.reasons.some(input => input.checked), false);
      assert.equal(app.run('events.size'), 0);
      assert.match(app.elements.counts.textContent, /2 unreviewed cases omitted/);
      assert.deepEqual(app.calls.map(call => call.url), ['cases.json']);
      assert.deepEqual(app.images, [row().raw, row().mask]);
      assert.equal(app.elements['raw-image'].src, row().raw);
      assert.equal(app.elements['raw-image'].style.width, '30px');
      assert.equal(app.elements['raw-image'].style.height, '20px');
      assert.equal(app.elements['mask-overlay'].style.maskImage, `url("${row().mask}")`);
      assert.equal(app.elements['mask-overlay'].style.opacity, .35);
      app.elements.zoom.value = '4'; await app.elements.zoom.emit('input');
      assert.equal(app.elements['raw-image'].style.width, '120px');
      assert.equal(app.elements['image-stage'].style.height, '80px');
      await app.document.emit('keydown', {key: 'm', target: {tagName: 'DIV'}, preventDefault() {}});
      assert.equal(app.elements['mask-overlay'].hidden, true);
      await app.document.emit('keydown', {key: 'ArrowRight', target: {tagName: 'TEXTAREA'}});
      assert.equal(app.elements['case-title'].textContent, 'case_01');
      await app.elements.next.click(); await tick();
      assert.equal(app.elements['case-title'].textContent, 'case_02');
      assert.equal(app.decisions.some(input => input.checked), false);
    """)


@pytest.mark.parametrize(
    "change",
    [
        "rows[0].specimen_id = 'hidden';",
        "rows[0].old_iod = 45;",
        "rows[0].raw = '../private/raw.png';",
        "rows[0].raw = 'https://example.test/raw.png';",
        "rows[0].mask = '/private/mask.png';",
        "rows[0].mask = 'images/%2e%2e/mask.png';",
        "rows[0].width = 1.5;",
        "rows[0].state_sha256 = 'unknown';",
        "rows.push({...rows[0]});",
    ],
)
def test_rejects_extra_metadata_invalid_state_or_unconfined_images(change):
    run_page(
        "const rows = [row()]; " + change + r"""
        const app = page(rows); await tick();
        assert.match(app.elements.message.textContent, /queue could not load/);
        assert.equal(app.images.length, 0);
        assert.equal(app.downloads.length, 0);
        assert.equal(app.elements.message.textContent.includes('hidden'), false);
        """
    )


def test_export_requires_explicit_reviewer_decision_and_reason():
    run_page(r"""
      const app = page(); await tick();
      await app.elements.export.click();
      assert.match(app.elements.message.textContent, /No authored decisions/);
      await app.choose('keep', ''); await app.elements.export.click();
      assert.match(app.elements.message.textContent, /at least one reason/);
      await app.choose('keep', 'complete_isolated_nucleus'); await app.elements.export.click();
      assert.match(app.elements.message.textContent, /reviewer name/);
      assert.equal(app.downloads.length, 0);
      app.elements.reviewer.value = 'Reviewer A';
      await app.choose('problem', 'other'); await app.elements.export.click();
      assert.match(app.elements.message.textContent, /Explain the other reason/);
      assert.equal(app.downloads.length, 0);
    """)


def test_export_is_exact_authored_event_and_immutable_in_session():
    run_page(r"""
      const app = page(); await tick();
      app.elements.reviewer.value = ' Reviewer A ';
      await app.choose('problem', 'incomplete_boundary', 'Faint edge <script> is excluded.');
      await app.elements.export.click();
      const payload = JSON.parse(await app.downloads[0].text());
      assert.equal(payload.length, 1);
      assert.deepEqual(Object.keys(payload[0]).sort(), ['case_id','state_sha256','reviewer','time','decision','reasons','note'].sort());
      assert.equal(payload[0].case_id, 'case_01');
      assert.equal(payload[0].state_sha256, row().state_sha256);
      assert.equal(payload[0].reviewer, 'Reviewer A');
      assert.match(payload[0].time, /^\d{4}-\d\d-\d\dT.*Z$/);
      assert.equal(payload[0].decision, 'problem');
      assert.deepEqual(payload[0].reasons, ['incomplete_boundary']);
      assert.equal(payload[0].note, 'Faint edge <script> is excluded.');
      assert.match(app.elements.message.textContent, /Omitted 1 unreviewed/);
      assert.equal(app.decisions.every(input => input.disabled), true);
      assert.equal(app.elements.reviewer.disabled, true);
      assert.equal(app.run('Object.isFrozen(events.get("case_01"))'), true);
      assert.equal(app.run('Object.isFrozen(events.get("case_01").reasons)'), true);
      await app.choose('keep', 'complete_isolated_nucleus', 'later draft');
      await app.elements.export.click();
      assert.deepEqual(JSON.parse(await app.downloads[1].text()), payload);
    """)


def test_incomplete_second_case_blocks_entire_export_without_partial_events():
    run_page(r"""
      const app = page(); await tick(); app.elements.reviewer.value = 'Reviewer A';
      await app.choose('keep', 'complete_isolated_nucleus');
      await app.elements.next.click(); await tick();
      await app.choose('unsure', ''); await app.elements.export.click();
      assert.equal(app.run('events.size'), 0);
      assert.equal(app.downloads.length, 0);
      assert.match(app.elements.message.textContent, /at least one reason/);
      await app.choose('unsure', 'uncertain_boundary'); await app.elements.export.click();
      assert.equal(JSON.parse(await app.downloads[0].text()).length, 2);
    """)


def test_cache_is_bound_to_packet_and_state_and_restores_only_drafts():
    run_page(r"""
      const storage = new Map(); const app = page(undefined, storage); await tick();
      app.elements.reviewer.value = 'Reviewer A';
      await app.choose('keep', 'complete_isolated_nucleus'); await app.elements.export.click();
      const keys = [...storage.keys()];
      assert.match(keys[0], /^nucleus-review-draft-v1:[a-f0-9]{64}:case_01:[a-f0-9]{64}$/);
      assert.equal(Object.hasOwn(JSON.parse(storage.get(keys[0])), 'time'), false);
      const restored = page(undefined, storage); await tick();
      assert.equal(restored.decisions[0].checked, true);
      assert.equal(restored.run('events.size'), 0);
      assert.equal(restored.decisions[0].disabled, false);
      assert.equal(restored.downloads.length, 0);
      const changed = page([row('case_01', 'b'.repeat(64)), row('case_02')], storage); await tick();
      assert.equal(changed.decisions.some(input => input.checked), false);
      const otherPacket = page([row(), row('case_03')], storage); await tick();
      assert.equal(otherPacket.decisions.some(input => input.checked), false);
    """)


@pytest.mark.parametrize("options", ["{badDimensions:true}", "{imageFailure:true}"])
def test_image_failure_prevents_review_and_export(options):
    run_page("const app = page(undefined, undefined, " + options + r"""); await tick();
      assert.equal(app.decisions.every(input => input.disabled), true);
      assert.equal(app.elements['image-stage'].hidden, true);
      assert.match(app.elements.message.textContent, /Image inspection failed/);
      app.elements.reviewer.value = 'Reviewer A';
      await app.choose('keep', 'complete_isolated_nucleus'); await app.elements.export.click();
      assert.equal(app.downloads.length, 0);
      assert.equal(app.run('events.size'), 0);
    """)


def test_storage_unavailable_retains_functional_explicit_export():
    run_page(r"""
      const app = page(undefined, undefined, {storageBlocked:true}); await tick();
      assert.match(app.elements['cache-status'].textContent, /storage is unavailable/);
      app.elements.reviewer.value = 'Reviewer A';
      await app.choose('unsure', 'uncertain_boundary'); await app.elements.export.click();
      assert.equal(JSON.parse(await app.downloads[0].text())[0].decision, 'unsure');
    """)


def test_markup_has_blank_decisions_native_mask_and_no_unblinding_or_editor():
    class Controls(HTMLParser):
        inputs: list[dict[str, str | None]]

        def __init__(self):
            super().__init__()
            self.inputs = []

        def handle_starttag(self, tag, attrs):
            if tag == "input":
                self.inputs.append(dict(attrs))

    source = PAGE.read_text()
    controls = Controls()
    controls.feed(source)
    decisions = [item for item in controls.inputs if item.get("name") == "decision"]
    assert {item["value"] for item in decisions} == {"keep", "problem", "unsure"}
    assert all("checked" not in item for item in decisions)
    assert all("checked" not in item for item in controls.inputs if item.get("name") == "reason")
    assert "mask-mode: luminance" in source
    assert "image-rendering: pixelated" in source
    for private_or_editor in [
        "specimen_id", "cohort", "process_id", "old_iod", "expected_dna", "source_path",
        "getImageData", "putImageData", "<canvas", "contenteditable=", "innerHTML",
    ]:
        assert private_or_editor not in source
