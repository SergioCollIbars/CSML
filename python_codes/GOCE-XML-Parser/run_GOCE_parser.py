import os
import shutil
import subprocess
from pathlib import Path

# === CONFIGURATION ===
script_dir = Path(__file__).resolve().parent
input_folder = script_dir / "input_files"
output_folder = script_dir / "output_files"
perl_script = script_dir / "hpf_eef_transform.pl"
keep_extensions = [".txt", ".csv", ".out"]  # adjust to match your output files

# === SETUP ===
output_folder.mkdir(exist_ok=True)
original_files = set(script_dir.glob("*"))  # Files in script directory before running

# === PROCESS EACH HDR + DBL PAIR ===
for hdr_file in input_folder.glob("*.HDR"):
    base_name = hdr_file.stem
    dbl_file = input_folder / (base_name + ".DBL")

    # Ensure both files exist
    if not dbl_file.exists():
        print(f"⚠️  Skipping {hdr_file.name} — missing DBL file.")
        continue

    print(f"🔄 Processing pair: {hdr_file.name}, {dbl_file.name}")

    # Record all files before run (recursively)
    before_files = {f.resolve() for f in script_dir.rglob("*")}

    # Run the Perl script with both files
    try:
        subprocess.run(["perl", str(perl_script), str(hdr_file), str(dbl_file)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error while processing {hdr_file.name}: {e}")
        continue

    # Record all files after run (recursively)
    after_files = {f.resolve() for f in script_dir.rglob("*")}
    new_files = after_files - before_files

    # Move new files to output folder
    for f in new_files:
        if f.is_file():
            destination = output_folder / f.name
            try:
                shutil.move(str(f), str(destination))
                print(f"📦 Moved: {f.name}")
            except Exception as e:
                print(f"⚠️  Failed to move {f.name}: {e}")

    print(f"✅ Finished {base_name}\n")

print("🏁 All pairs processed.")
