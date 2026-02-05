# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['Matrix_editor.pyw'],
    pathex=[],
    binaries=[('D:\\pymol\\Library\\bin\\tcl86t.dll', '.'), ('D:\\pymol\\Library\\bin\\tk86t.dll', '.')],
    datas=[('D:\\pymol\\Library\\lib\\tcl8.6', 'tcl\\tcl8.6'), ('D:\\pymol\\Library\\lib\\tk8.6', 'tcl\\tk8.6')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='MatrixEditor_debug',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
