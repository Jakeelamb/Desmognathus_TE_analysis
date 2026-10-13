/* Figure review notes: project saves with revision checks and local draft backups. */
(() => {
  "use strict";

  function start() {
    const summary = document.querySelector("#notes-summary");
    const saveButton = document.querySelector("#save-notes");
    const exportButton = document.querySelector("#export-notes");
    if (!summary || !saveButton || !exportButton) return;

    const prefix = `desmognathus.figure-review.notes.v1:${location.pathname}:`;
    const states = [];
    let removedNotes = [];
    let ready = false;
    let loading = false;
    let loadFailed = false;

    function entry(value) {
      if (!value || typeof value.text !== "string" || value.text.length > 50000 ||
          !Number.isSafeInteger(value.revision) || value.revision < 0) {
        throw new Error("Invalid notes response");
      }
      return value;
    }

    function backup(state) {
      try {
        if (state.dirty || state.saving || state.conflict) {
          localStorage.setItem(prefix + state.id, JSON.stringify({
            text: state.input.value,
            title: state.title,
            base_revision: state.baseRevision,
          }));
        } else {
          localStorage.removeItem(prefix + state.id);
        }
        state.backupOk = true;
      } catch (_) {
        state.backupOk = false;
      }
    }

    function render() {
      const removedCount = removedNotes.filter(note => note.text.trim()).length;
      let count = removedCount;
      let pending = 0;
      let conflicts = 0;
      let failures = 0;
      for (const state of states) {
        if (state.input.value.trim()) count += 1;
        if (state.dirty || state.saving) pending += 1;
        if (state.conflict) conflicts += 1;
        if (state.error) failures += 1;
        const local = state.backupOk
          ? "Draft saved in this browser."
          : "Browser backup unavailable; keep this page open or export your notes.";
        let message;
        if (state.conflict) {
          message = `A different version is saved in the project. ${local}`;
        } else if (state.saving) {
          message = `Saving to project… ${local}`;
        } else if (state.dirty) {
          message = `${state.error || !ready ? "Project save unavailable. " : ""}${local}`;
        } else if (!ready) {
          message = loading ? "Loading project notes…" : "Project notes unavailable. You can write a browser draft.";
        } else {
          message = state.saved.revision > 0 ? "Saved to project" : "No notes yet";
        }
        state.status.textContent = message;
        state.conflictBox.hidden = !state.conflict;
        state.remote.textContent = state.conflict
          ? state.saved.text || "(The project note is empty.)" : "";
      }
      const pieces = [`${count} ${count === 1 ? "figure/panel has" : "figures/panels have"} notes`];
      if (removedCount) pieces.push(`${removedCount} ${removedCount === 1 ? "note on a removed panel" : "notes on removed panels"}`);
      if (pending) pieces.push(`${pending} unsaved`);
      if (conflicts) pieces.push(`${conflicts} ${conflicts === 1 ? "conflict" : "conflicts"} to resolve`);
      if (loading) pieces.push("loading project notes…");
      else if (loadFailed || failures) pieces.push("project save unavailable; use Save now to retry");
      else if (ready && !pending) pieces.push("all notes saved to project");
      summary.textContent = pieces.join(" · ");
      saveButton.disabled = loading || (!pending && ready);
      exportButton.disabled = count === 0;
    }

    async function request(url, options = {}) {
      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), 8000);
      try {
        const response = await fetch(url, { ...options, signal: controller.signal, cache: "no-store" });
        const body = await response.json();
        if (!response.ok && response.status !== 409) throw new Error("Project request failed");
        return { response, body };
      } finally {
        clearTimeout(timeout);
      }
    }

    function queue(state) {
      clearTimeout(state.timer);
      if (ready && state.dirty && !state.conflict) {
        state.timer = setTimeout(() => save(state), 600);
      }
    }

    function save(state) {
      clearTimeout(state.timer);
      if (state.flight) return state.flight;
      if (!ready || !state.dirty || state.conflict) return Promise.resolve();
      state.error = false;
      state.flight = (async () => {
        try {
          // Drain edits made during a request only after its revision is acknowledged.
          while (state.dirty && !state.conflict) {
            const text = state.input.value;
            state.saving = true;
            backup(state);
            render();
            const { response, body } = await request(`/api/notes/${encodeURIComponent(state.id)}`, {
              method: "PUT",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ text, title: state.title, base_revision: state.baseRevision }),
            });
            const saved = entry(response.status === 409 && body.entry === null
              ? { text: "", revision: 0 } : body.entry);
            if (response.status !== 409 && saved.text !== text) {
              throw new Error("Project returned a different note");
            }
            state.saved = saved;
            state.uncertain = false;
            state.saving = false;
            state.dirty = state.input.value !== saved.text;
            state.conflict = response.status === 409 && state.dirty;
            if (!state.conflict) state.baseRevision = saved.revision;
            backup(state);
          }
        } catch (_) {
          state.error = true;
          // A lost acknowledgement may still mean the server committed the write.
          state.uncertain = true;
          state.dirty = true;
        } finally {
          state.saving = false;
          state.flight = null;
          backup(state);
          render();
        }
      })();
      return state.flight;
    }

    async function load() {
      if (loading) return;
      loading = true;
      render();
      try {
        const { response, body } = await request("/api/notes");
        if (!response.ok || body.schema_version !== 1 || !body.notes ||
            typeof body.notes !== "object" || Array.isArray(body.notes)) {
          throw new Error("Invalid notes response");
        }
        // Validate the entire response before applying any entries.
        const savedById = new Map(Object.entries(body.notes).map(([id, value]) => [id, entry(value)]));
        const saved = states.map(state => savedById.get(state.id) || { text: "", revision: 0 });
        const currentIds = new Set(states.map(state => state.id));
        removedNotes = [...savedById].filter(([id]) => !currentIds.has(id))
          .map(([id, note]) => ({ ...note, id }));
        states.forEach((state, index) => {
          state.saved = saved[index];
          state.error = false;
          state.uncertain = false;
          if (!state.dirty) state.input.value = state.saved.text;
          state.dirty = state.input.value !== state.saved.text;
          state.conflict = state.dirty && state.baseRevision !== state.saved.revision;
          if (!state.conflict) state.baseRevision = state.saved.revision;
          backup(state);
        });
        ready = true;
        loadFailed = false;
      } catch (_) {
        loadFailed = true;
      } finally {
        loading = false;
        for (const state of states) {
          state.input.disabled = false;
          queue(state);
        }
        render();
      }
    }

    for (const input of document.querySelectorAll(".note-input[data-note-id]")) {
      const section = input.closest("section[id]");
      const state = {
        id: input.dataset.noteId,
        title: section.querySelector("h2").textContent.trim(),
        input,
        status: section.querySelector(".note-status"),
        conflictBox: section.querySelector(".note-conflict"),
        remote: section.querySelector(".remote-note"),
        saved: { text: "", revision: 0 },
        baseRevision: 0,
        dirty: false,
        saving: false,
        conflict: false,
        error: false,
        uncertain: false,
        backupOk: true,
        flight: null,
        timer: null,
      };
      try {
        const draft = JSON.parse(localStorage.getItem(prefix + state.id));
        if (draft && typeof draft.text === "string" && draft.text.length <= 50000 &&
            Number.isSafeInteger(draft.base_revision) && draft.base_revision >= 0) {
          input.value = draft.text;
          state.baseRevision = draft.base_revision;
          state.dirty = true;
        }
      } catch (_) {
        state.backupOk = false;
      }
      states.push(state);
      input.addEventListener("input", () => {
        state.dirty = state.saving || state.uncertain || input.value !== state.saved.text;
        if (!state.dirty) {
          state.conflict = false;
          state.baseRevision = state.saved.revision;
        }
        state.error = false;
        backup(state);
        render();
        queue(state);
      });
      input.addEventListener("blur", () => save(state));
      section.querySelector('[data-note-action="focus"]').addEventListener("click", () => input.focus());
      section.querySelector('[data-note-action="use-saved"]').addEventListener("click", () => {
        input.value = state.saved.text;
        state.baseRevision = state.saved.revision;
        state.conflict = false;
        state.dirty = false;
        state.error = false;
        backup(state);
        render();
      });
      section.querySelector('[data-note-action="keep-draft"]').addEventListener("click", () => {
        state.baseRevision = state.saved.revision;
        state.conflict = false;
        state.dirty = input.value !== state.saved.text;
        backup(state);
        save(state);
        render();
      });
    }

    saveButton.addEventListener("click", async () => {
      if (!ready) await load();
      await Promise.all(states.map(save));
    });
    exportButton.addEventListener("click", () => {
      const lines = ["# Figure review notes", "", `Exported ${new Date().toISOString()}`, ""];
      for (const state of states) {
        if (!state.input.value.trim()) continue;
        const url = new URL(location.href);
        url.hash = state.id;
        lines.push(`## ${state.title}`, "", `[View figure or panel](${url.href})`, "");
        if (state.dirty || state.saving || state.conflict) {
          lines.push("*This note includes a draft that has not been saved to the project.*", "");
        }
        lines.push(state.input.value, "");
      }
      for (const note of removedNotes) {
        if (!note.text.trim()) continue;
        const title = typeof note.title === "string" && note.title.trim() ? note.title : note.id;
        lines.push(`## ${title}`, "", "*Saved note. Panel removed from current gallery.*", "", note.text, "");
      }
      const url = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/markdown;charset=utf-8" }));
      const link = document.createElement("a");
      link.href = url;
      link.download = `figure-review-notes-${new Date().toISOString().slice(0, 10)}.md`;
      document.body.appendChild(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    });
    window.addEventListener("beforeunload", event => {
      if (states.some(state => state.dirty || state.saving || state.conflict)) {
        event.preventDefault();
        event.returnValue = "";
      }
    });
    load();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start, { once: true });
  else start();
})();
