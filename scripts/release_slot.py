from filelock import FileLock
import sys

if len(sys.argv) != 2:
    raise ValueError("Usage: release_slot.py <slot>")

slot = sys.argv[1]
SLOTS_FILE = "gpu_slots.txt"
LOCK_FILE = "gpu_slots.lock"

with FileLock(LOCK_FILE):
    with open(SLOTS_FILE, "r") as f:
        slots = f.read().split()

    slots.append(slot)

    with open(SLOTS_FILE, "w") as f:
        f.write("\n".join(slots))
