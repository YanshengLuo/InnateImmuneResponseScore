#!/usr/bin/env python3
"""
Design table editor (TSV/CSV) with a simple UI.

What it does:
- Pick a .tsv/.csv file (e.g., *_design.tsv)
- Load into a table view
- Click a cell to edit it (or use controls to set an entire column or conditional replace)
- Save as TSV/CSV

Works on Windows/macOS/Linux with standard Python + tkinter.
"""

import os
import sys
import csv
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# ---------- Utilities ----------

def detect_delimiter(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    if ext == ".tsv":
        return "\t"
    if ext == ".csv":
        return ","
    # fallback: sniff
    with open(path, "r", newline="", encoding="utf-8") as f:
        sample = f.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except Exception:
        return "\t"

def read_table(path: str):
    delim = detect_delimiter(path)
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=delim)
        rows = list(reader)
    if not rows:
        return [], [], delim
    header = rows[0]
    data = rows[1:]
    # Normalize row widths
    max_cols = len(header)
    for r in data:
        if len(r) < max_cols:
            r.extend([""] * (max_cols - len(r)))
        elif len(r) > max_cols:
            # If file has extra columns, extend header
            extra = len(r) - max_cols
            header.extend([f"extra_{i+1}" for i in range(extra)])
            max_cols = len(header)
    return header, data, delim

def write_table(path: str, header, data, delim: str):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=delim)
        writer.writerow(header)
        writer.writerows(data)

def safe_int(s: str):
    try:
        return int(s)
    except Exception:
        return None

def looks_numeric_col(values):
    # if most non-empty values are numeric-like
    non_empty = [v for v in values if str(v).strip() != ""]
    if not non_empty:
        return False
    numeric = 0
    for v in non_empty:
        try:
            float(str(v))
            numeric += 1
        except Exception:
            pass
    return numeric / len(non_empty) >= 0.8

def coerce_value(col_is_numeric: bool, new_val: str):
    if not col_is_numeric:
        return new_val
    try:
        # keep integers pretty if possible
        fv = float(new_val)
        if fv.is_integer():
            return str(int(fv))
        return str(fv)
    except Exception:
        # if numeric column but value not numeric, keep as-is (user might want strings)
        return new_val

# ---------- UI App ----------

class TableEditorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("IMRS Design Editor (TSV/CSV)")
        self.geometry("1100x700")

        self.file_path = None
        self.delim = "\t"
        self.header = []
        self.data = []   # list[list[str]]

        self._build_ui()

    def _build_ui(self):
        # Top bar
        top = ttk.Frame(self, padding=10)
        top.pack(fill="x")

        ttk.Button(top, text="Open File…", command=self.open_file).pack(side="left")
        ttk.Button(top, text="Save", command=self.save_file).pack(side="left", padx=(8,0))
        ttk.Button(top, text="Save As…", command=self.save_file_as).pack(side="left", padx=(8,0))

        self.path_label = ttk.Label(top, text="No file loaded", foreground="#444")
        self.path_label.pack(side="left", padx=12)

        # Main split
        main = ttk.PanedWindow(self, orient="horizontal")
        main.pack(fill="both", expand=True, padx=10, pady=(0,10))

        left = ttk.Frame(main)
        right = ttk.Frame(main, padding=10)
        main.add(left, weight=3)
        main.add(right, weight=1)

        # Table (Treeview)
        self.tree = ttk.Treeview(left, show="headings")
        vsb = ttk.Scrollbar(left, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(left, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.pack(fill="both", expand=True, side="top")
        vsb.pack(fill="y", side="right")
        hsb.pack(fill="x", side="bottom")

        self.tree.bind("<Double-1>", self.on_double_click)

        # Right controls
        ttk.Label(right, text="Cell edit", font=("Segoe UI", 11, "bold")).pack(anchor="w")

        self.selected_info = ttk.Label(right, text="Select a cell (double-click) to edit.", wraplength=320)
        self.selected_info.pack(anchor="w", pady=(6,10))

        cell_box = ttk.Frame(right)
        cell_box.pack(fill="x")

        ttk.Label(cell_box, text="Row (1-based):").grid(row=0, column=0, sticky="w")
        self.row_entry = ttk.Entry(cell_box, width=10)
        self.row_entry.grid(row=0, column=1, sticky="w", padx=(8,0))

        ttk.Label(cell_box, text="Column:").grid(row=1, column=0, sticky="w", pady=(6,0))
        self.col_combo = ttk.Combobox(cell_box, state="readonly", values=[])
        self.col_combo.grid(row=1, column=1, sticky="we", padx=(8,0), pady=(6,0))

        ttk.Label(cell_box, text="New value:").grid(row=2, column=0, sticky="w", pady=(6,0))
        self.val_entry = ttk.Entry(cell_box)
        self.val_entry.grid(row=2, column=1, sticky="we", padx=(8,0), pady=(6,0))

        cell_box.columnconfigure(1, weight=1)

        ttk.Button(right, text="Apply to this cell", command=self.apply_cell).pack(anchor="w", pady=(8,14))

        ttk.Separator(right).pack(fill="x", pady=10)

        # Column set
        ttk.Label(right, text="Set entire column", font=("Segoe UI", 11, "bold")).pack(anchor="w")
        colset = ttk.Frame(right)
        colset.pack(fill="x", pady=(6,0))

        ttk.Label(colset, text="Column:").grid(row=0, column=0, sticky="w")
        self.set_col_combo = ttk.Combobox(colset, state="readonly", values=[])
        self.set_col_combo.grid(row=0, column=1, sticky="we", padx=(8,0))

        ttk.Label(colset, text="Value:").grid(row=1, column=0, sticky="w", pady=(6,0))
        self.set_col_val = ttk.Entry(colset)
        self.set_col_val.grid(row=1, column=1, sticky="we", padx=(8,0), pady=(6,0))

        colset.columnconfigure(1, weight=1)

        ttk.Button(right, text="Set column", command=self.set_column).pack(anchor="w", pady=(8,14))

        ttk.Separator(right).pack(fill="x", pady=10)

        # Conditional replace (simple)
        ttk.Label(right, text="Conditional replace", font=("Segoe UI", 11, "bold")).pack(anchor="w")
        cond = ttk.Frame(right)
        cond.pack(fill="x", pady=(6,0))

        ttk.Label(cond, text="Where col:").grid(row=0, column=0, sticky="w")
        self.where_col = ttk.Combobox(cond, state="readonly", values=[])
        self.where_col.grid(row=0, column=1, sticky="we", padx=(8,0))

        ttk.Label(cond, text="equals:").grid(row=1, column=0, sticky="w", pady=(6,0))
        self.where_val = ttk.Entry(cond)
        self.where_val.grid(row=1, column=1, sticky="we", padx=(8,0), pady=(6,0))

        ttk.Label(cond, text="Set col:").grid(row=2, column=0, sticky="w", pady=(6,0))
        self.to_col = ttk.Combobox(cond, state="readonly", values=[])
        self.to_col.grid(row=2, column=1, sticky="we", padx=(8,0), pady=(6,0))

        ttk.Label(cond, text="to value:").grid(row=3, column=0, sticky="w", pady=(6,0))
        self.to_val = ttk.Entry(cond)
        self.to_val.grid(row=3, column=1, sticky="we", padx=(8,0), pady=(6,0))

        cond.columnconfigure(1, weight=1)

        ttk.Button(right, text="Apply conditional replace", command=self.conditional_replace).pack(anchor="w", pady=(8,0))

        # Footer hint
        ttk.Label(
            right,
            text="Tip: Save As if you want to keep the original.\nDouble-click a cell to prefill row/col.",
            foreground="#555",
            wraplength=320
        ).pack(anchor="w", pady=(14,0))

    # ---------- File actions ----------

    def open_file(self):
        path = filedialog.askopenfilename(
            title="Select design TSV/CSV",
            filetypes=[("TSV files", "*.tsv"), ("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            header, data, delim = read_table(path)
        except Exception as e:
            messagebox.showerror("Open failed", str(e))
            return

        self.file_path = path
        self.delim = delim
        self.header = header
        self.data = data

        self.path_label.configure(text=os.path.basename(path))
        self._load_into_tree()
        self._refresh_column_controls()

    def save_file(self):
        if not self.file_path:
            messagebox.showinfo("No file", "Open a file first.")
            return
        try:
            write_table(self.file_path, self.header, self.data, self.delim)
            messagebox.showinfo("Saved", f"Saved:\n{self.file_path}")
        except Exception as e:
            messagebox.showerror("Save failed", str(e))

    def save_file_as(self):
        if not self.file_path:
            messagebox.showinfo("No file", "Open a file first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save As",
            initialfile=os.path.basename(self.file_path),
            defaultextension=os.path.splitext(self.file_path)[1],
            filetypes=[("TSV files", "*.tsv"), ("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not path:
            return
        delim = "\t" if path.lower().endswith(".tsv") else ("," if path.lower().endswith(".csv") else self.delim)
        try:
            write_table(path, self.header, self.data, delim)
            messagebox.showinfo("Saved", f"Saved:\n{path}")
        except Exception as e:
            messagebox.showerror("Save failed", str(e))

    # ---------- Table rendering ----------

    def _load_into_tree(self):
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = self.header

        # Configure headings + widths
        for col in self.header:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=max(90, min(240, 12 * len(col))), stretch=True)

        # Insert rows (prefix with row number in iid for stable access)
        for i, row in enumerate(self.data, start=1):
            # Ensure correct length
            if len(row) < len(self.header):
                row = row + [""] * (len(self.header) - len(row))
                self.data[i-1] = row
            self.tree.insert("", "end", iid=str(i), values=row)

    def _refresh_column_controls(self):
        cols = self.header
        self.col_combo["values"] = cols
        self.set_col_combo["values"] = cols
        self.where_col["values"] = cols
        self.to_col["values"] = cols

        # default selections
        if cols:
            self.col_combo.set(cols[0])
            self.set_col_combo.set(cols[0])
            self.where_col.set(cols[0])
            self.to_col.set(cols[0])

    # ---------- Cell editing ----------

    def on_double_click(self, event):
        if not self.header:
            return
        region = self.tree.identify("region", event.x, event.y)
        if region != "cell":
            return

        row_id = self.tree.identify_row(event.y)
        col_id = self.tree.identify_column(event.x)  # like '#3'
        if not row_id or not col_id:
            return

        r = int(row_id)  # 1-based
        c_index = int(col_id.replace("#", "")) - 1
        if c_index < 0 or c_index >= len(self.header):
            return

        col_name = self.header[c_index]
        current_val = self.data[r-1][c_index]

        self.row_entry.delete(0, tk.END)
        self.row_entry.insert(0, str(r))
        self.col_combo.set(col_name)
        self.val_entry.delete(0, tk.END)
        self.val_entry.insert(0, current_val)

        self.selected_info.configure(text=f"Selected: row {r}, column '{col_name}'\nCurrent value: {current_val}")

    def apply_cell(self):
        if not self.header:
            messagebox.showinfo("No file", "Open a file first.")
            return

        r = safe_int(self.row_entry.get().strip())
        if r is None or r < 1 or r > len(self.data):
            messagebox.showerror("Invalid row", f"Row must be between 1 and {len(self.data)}")
            return

        col = self.col_combo.get().strip()
        if col not in self.header:
            messagebox.showerror("Invalid column", "Choose a valid column.")
            return
        c = self.header.index(col)

        new_val_raw = self.val_entry.get()
        col_vals = [row[c] for row in self.data]
        col_is_num = looks_numeric_col(col_vals)

        new_val = coerce_value(col_is_num, new_val_raw)

        self.data[r-1][c] = new_val
        self.tree.item(str(r), values=self.data[r-1])

        self.selected_info.configure(text=f"Updated: row {r}, column '{col}' -> {new_val}")

    # ---------- Bulk operations ----------

    def set_column(self):
        if not self.header:
            messagebox.showinfo("No file", "Open a file first.")
            return

        col = self.set_col_combo.get().strip()
        if col not in self.header:
            messagebox.showerror("Invalid column", "Choose a valid column.")
            return
        c = self.header.index(col)
        val_raw = self.set_col_val.get()

        col_vals = [row[c] for row in self.data]
        col_is_num = looks_numeric_col(col_vals)
        val = coerce_value(col_is_num, val_raw)

        for i in range(len(self.data)):
            self.data[i][c] = val
            self.tree.item(str(i+1), values=self.data[i])

        messagebox.showinfo("Done", f"Set column '{col}' to '{val}' for {len(self.data)} rows.")

    def conditional_replace(self):
        if not self.header:
            messagebox.showinfo("No file", "Open a file first.")
            return

        where_col = self.where_col.get().strip()
        to_col = self.to_col.get().strip()
        if where_col not in self.header or to_col not in self.header:
            messagebox.showerror("Invalid column", "Choose valid columns.")
            return

        where_val = self.where_val.get()
        to_val_raw = self.to_val.get()

        wc = self.header.index(where_col)
        tc = self.header.index(to_col)

        # infer numeric type for target col
        target_vals = [row[tc] for row in self.data]
        target_is_num = looks_numeric_col(target_vals)
        to_val = coerce_value(target_is_num, to_val_raw)

        hit = 0
        for i, row in enumerate(self.data):
            if str(row[wc]) == str(where_val):
                row[tc] = to_val
                self.tree.item(str(i+1), values=row)
                hit += 1

        messagebox.showinfo("Done", f"Replaced {hit} rows where {where_col} == '{where_val}' : set {to_col} = '{to_val}'.")

# ---------- Run ----------

if __name__ == "__main__":
    app = TableEditorApp()
    app.mainloop()
