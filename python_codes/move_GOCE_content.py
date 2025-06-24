import os
import shutil

def move_GO_CONS_contents(root_dir, destination_dir, delete_empty_folders=False):
    """
    Move contents of all folders with 'GO_CONS' in the name to a target directory.

    Parameters:
    - root_dir: path to search for GO_CONS folders (searches recursively)
    - destination_dir: where to move the contents
    - delete_empty_folders: if True, deletes the source folder after moving
    """
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    for dirpath, dirnames, _ in os.walk(root_dir):
        for dirname in dirnames:
            if "GO_CONS" in dirname:
                source_folder = os.path.join(dirpath, dirname)
                print(f"Found: {source_folder}")

                for item in os.listdir(source_folder):
                    item_path = os.path.join(source_folder, item)
                    if os.path.isfile(item_path):
                        shutil.move(item_path, os.path.join(destination_dir, item))
                        print(f"Moved {item} → {destination_dir}")

                if delete_empty_folders and not os.listdir(source_folder):
                    os.rmdir(source_folder)
                    print(f"Deleted empty folder: {source_folder}")

# Example usage
root_directory = '/Users/sergiocollibars/Downloads'
destination_directory = '/Users/sergiocollibars/Desktop/CSML/python_codes/GOCE-XML-Parser/input_files'

move_GO_CONS_contents(root_directory, destination_directory, delete_empty_folders=True)
