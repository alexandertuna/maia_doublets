# Run directory

This is how I run the code:

```bash
# Go to the run/ directory
cd run/

# Create a virtual environment (once)
python -m venv env_ubuntu24

# Activate the environment (on each login)
source env_ubuntu24/bin/activate

# pip install the repo
pip install -e ..

# Run!
maia_doublets -h
```
