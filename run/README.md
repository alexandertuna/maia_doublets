# Run directory

This is how I run the code:

```bash
# Go to the run/ directory
cd run/

# Create a virtual environment (once)
python -m venv env_ubuntu24

# Set up the MuColl stack and activate the environment (on each login)
source env.sh

# pip install the repo
pip install -e ..

# Run!
maia_doublets -h
```
