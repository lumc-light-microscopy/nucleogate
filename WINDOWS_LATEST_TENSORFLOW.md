# Running the latest TensorFlow on Windows (via WSL2)

`environment.yml` pins TensorFlow 2.10 and Python 3.10 because **2.10 was the
last TensorFlow release with native-Windows GPU support**. If you want current
TensorFlow *and* a GPU on a Windows machine, the supported route is WSL2 —
a real Linux environment running on your Windows PC, with access to the same
GPU. This page is the short version.

**You probably don't need this if:**

- You're happy on the pinned TF 2.10 stack — it works, it's just frozen.
- You don't have an NVIDIA GPU. Current TensorFlow runs natively on Windows
  CPU-only: `pip install tensorflow` and no WSL2 needed.
- You'd rather not maintain a Linux install — use `environment-cellpose.yml`
  instead and run `--model cellpose`. PyTorch still ships native-Windows GPU
  builds, so you get a current stack with none of this setup.

Only the segmentation backend cares. Everything else in `omero_nuclei_flow.py`
runs identically either way.

---

## 1. Prerequisites

- Windows 11, or Windows 10 build 19044+ (21H2, the Nov 2021 update).
- An NVIDIA GPU with a current driver **installed on Windows**.
- **Do not install an NVIDIA driver inside WSL.** WSL2 uses the Windows
  driver. Installing a Linux driver in the guest is the single most common way
  people break GPU passthrough.

## 2. Install WSL2

In PowerShell **as administrator**:

```powershell
wsl --install -d Ubuntu-24.04
```

Reboot, then open "Ubuntu" from the Start menu. It'll ask you to create a Linux
username and password (unrelated to your Windows login). Everything from here
runs inside that Ubuntu terminal.

Check the GPU is visible:

```bash
nvidia-smi
```

If that lists your card, passthrough works. If it doesn't, stop and fix that
before installing anything else.

## 3. Install Miniforge inside WSL

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

Accept the defaults, then close and reopen the Ubuntu terminal.

## 4. Create the environment

```bash
conda create -n omeroflow-wsl -c conda-forge python=3.11 \
    zeroc-ice=3.6.5 omero-py "numpy<2" pandas scipy scikit-image matplotlib \
    scikit-learn pyarrow openpyxl "setuptools<81"
conda activate omeroflow-wsl
```

Then TensorFlow with CUDA, from pip (this pulls the matching CUDA libraries
itself — you don't install a CUDA toolkit separately):

```bash
pip install "tensorflow[and-cuda]"
pip install stardist tf-keras
```

Check the current syntax at <https://www.tensorflow.org/install/pip> if that
fails; it's the one command in this guide that has changed over time.

## 5. Verify

```bash
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
python omero_nuclei_flow.py --check-env
```

You want a non-empty GPU list and `environment looks consistent`.

## 6. Run

```bash
cd /mnt/c/OmeroFlow          # your Windows folder, seen from Linux
python omero_nuclei_flow.py --host omero.myinstitute.org --user you \
    --screen-id 1234 --outdir /mnt/c/OmeroFlow/results
```

Windows drives appear under `/mnt/c`, `/mnt/d` and so on, so results written
there are immediately visible in Windows Explorer.

---

## Three things that catch people out

### Keras 3

**This is the main trap.** From TensorFlow 2.16 onward, `pip install
tensorflow` installs **Keras 3**, and StarDist only supports Keras 2. That's why
step 4 installs `tf-keras` (the legacy Keras 2 package). The script detects it
and sets `TF_USE_LEGACY_KERAS=1` automatically before importing TensorFlow, and
logs that it did so:

```
tf-keras found -> TF_USE_LEGACY_KERAS=1 (StarDist needs Keras 2)
```

If you skip `tf-keras`, `--check-env` will tell you. Symptoms otherwise are
obscure `AttributeError`s from inside StarDist about missing Keras attributes.

### numpy 2

StarDist is not numpy-2 clean. Step 4 pins `numpy<2` for this reason. If
something later upgrades numpy, expect `numpy.core.multiarray failed to import`
or similar — `--check-env` flags this too.

### Disk speed and VPNs

- Reading image data across `/mnt/c` is slow. That doesn't matter here, since
  images come over the network from OMERO, but if you ever process local files,
  copy them into the Linux filesystem (`~/data`) first.
- **If your OMERO server is behind an institutional VPN running on Windows,
  WSL2 may not see it.** On Windows 11 23H2+, create `C:\Users\<you>\.wslconfig`
  containing:

  ```ini
  [wsl2]
  networkingMode=mirrored
  ```

  then run `wsl --shutdown` in PowerShell and reopen Ubuntu. This makes WSL2
  share the Windows network stack, VPN included. Test with
  `ping omero.myinstitute.org` before blaming the script.

---

## Going back

Nothing here touches your Windows install. The pinned Windows environment
(`conda activate omeroflow` in the Miniforge Prompt) still works exactly as
before, and both can coexist. To remove WSL2 entirely:
`wsl --unregister Ubuntu-24.04` in PowerShell.
