import sys
from pathlib import Path
import csv
import os

if len(sys.argv) != 4:
    print("This program needs 3 arguments:")
    print("- Project name\n- Manifest file with three columns: sample-id, absolute-filepath, direction\n- Sample name")
    sys.exit()

script_path, project_name, manifest_path, sample_name = sys.argv

folders = [
    "results/" + project_name + "/samples/" + sample_name + "/rawdata/forward_reads",
    "results/" + project_name + "/samples/" + sample_name + "/rawdata/reverse_reads"
]

for folder in folders:
    Path(folder).mkdir(parents=True, exist_ok=True)
    print("Created folder %s" % (folder))

with open(manifest_path) as csv_file:
    reader = csv.reader(csv_file)
    next(reader)
    for row in reader:
        sample_id, file_path, direction = row

        if sample_id == sample_name:
            target_file = "fw.fastq.gz" if direction == "forward" else "rv.fastq.gz"
            target = os.path.abspath("results/" + project_name + "/samples/" + sample_id + "/rawdata/" + direction + "_reads/" + target_file)
            source = os.path.abspath(file_path)
            if Path(target).exists() or Path(target).is_symlink():
                Path(target).unlink()
                print("Removed existing symlink %s" % (target))
            os.symlink(source, target)
            print("Created symlink %s -> %s" % (target, source))
