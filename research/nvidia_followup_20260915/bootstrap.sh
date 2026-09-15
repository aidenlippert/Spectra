set -eu
mkdir -p /home/ubuntu/spectra-nvidia
cd /home/ubuntu/spectra-nvidia
python3 -m pip install --user uv
/home/ubuntu/.local/bin/uv python install 3.12
/home/ubuntu/.local/bin/uv venv --python 3.12 .venv
/home/ubuntu/.local/bin/uv pip install --python .venv/bin/python 'numpy==2.2.6' 'scipy==1.18.1' 'cupy-cuda12x[ctk]==14.2.0' 'cuquantum-python-cu12' 'nvmath-python[cu12]' 'nvidia-cudss-cu12' 'quimb==1.15.0' numba networkx threadpoolctl python-flint pytest
.venv/bin/python -m pip freeze > installed.txt 2>&1 || /home/ubuntu/.local/bin/uv pip freeze --python .venv/bin/python > installed.txt
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader > gpu.txt
sudo apt-get update -qq
sudo apt-get install -y -qq libsuitesparse-dev
/home/ubuntu/.local/bin/uv pip install --python .venv/bin/python 'scikit-sparse==0.4.16'
/home/ubuntu/.local/bin/uv pip freeze --python .venv/bin/python > installed.txt
