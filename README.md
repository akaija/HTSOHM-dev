#HTSOHM-dev

Master Branch Status: [![Build Status](https://travis-ci.org/WilmerLab/HTSOHM-dev.svg?branch=master)](https://travis-ci.org/WilmerLab/HTSOHM-dev)

This is a package which can be used to map the structure-property space   
spanned by hypothetical, porous <i>pseudo materials</i>. This package contains   
modules that manage randomly-generating libraries of pseudo-materials,   
distributing many simulations across a computing cluster, and selectively   
mutating rare materials--in an iterative process.

## Getting started
This package was built to be run locally, or on a computing cluster. To run on   
a cluster, the following must be installed:   
  - Torque/PBS (qsub)   
  - Python 3.5.1     
In any case, following Python packages are required:   
  - pyyaml   
  - psycopg2   
  - numpy   
  - SQLAlchemy   
  - RASPA2 (git+https://github.com/WilmerLab/raspa2)   
  
## First time setup
```
# If using environment modules    
module purge   
module load python/3.5.1/environment/module    

# Clone HTSOHM repo    
git clone https://github.com/WilmerLab/HTSOHM-dev.git    

# Create virual environment for Python packages    
mkdir -p ~/venv    
pyvenv ~/venv/htsohm    
source ~/venv/htsohm/bin/activate    

# Install Python packages   
cd path/to/htsohm/repo    
pip install -r requirements.txt    
```
You will need a database.yaml file for HTSOHM to use to store results. You    
can use the database.sample.yaml file (just copy it to database.yaml) if you    
just want to use a local SQLite database. Otherwise, enter in your    
connection_string into the database.yaml file per the format of the example    
file.    
  
You will also need configuration file to specify the run-parameters. You can    
use the htsohm.sample.yaml file (just copy it and change the values specified    
within).    
  
## Running HTSOHM
### Activating the htsohm virtual environment:
```
module purge
module load python/3.5.1/environment/module
source ~/venv/htsohm/bin/activate
cd path/to/htsohm/repo
```
### Starting a run:    
```
./hts.py start path/to/config      
```
### Launching a worker locally:    
```
./hts.py launch_worker run_id    
```
### Launching workers on a cluster:      
```
qsub -v='run_id' launch_workers_qsub.sh   
```
When running on a cluster, we recommend using GNU screen or Tmux; one screen    
running the htsohm virtual environment, and one to be used to access cluster    
resources such as PBS/Torque.

### Install to virtual environment:
```
cd /path/to/htsohm/repo
~/venv/htsohm/bin/python setup.py install
```

### Visualising pseudomaterials
```
blender empty_scene.blend --background --python blender_script.py <material uuid>
```
Note: Blender is a prerequisite for visualisation.  

  
Please send questions/comments/concerns to `ark111@pitt.edu`.
  
## License
  
HTSOHM-dev and related modules released under the MIT License 2016.    
  
## Development environment

```
3.5.1 |Anaconda 2.4.1 (64-bit)| (default, Dec  7 2015, 11:16:01)
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]
Bash 4.1.2(1)-release
RHEL 6.6	2.6.32-358
psql (PostgreSQL) 8.4.20
```
