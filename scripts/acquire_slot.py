from filelock import FileLock

SLOTS_FILE = "gpu_slots.txt"
LOCK_FILE = "gpu_slots.lock"

with FileLock(LOCK_FILE):
    with open(SLOTS_FILE, "r") as f:
        slots = f.read().split()

    if not slots:
        raise RuntimeError("No GPU slots available!")

    slot = slots.pop(0)

    with open(SLOTS_FILE, "w") as f:
        f.write("\n".join(slots))

print(slot)
