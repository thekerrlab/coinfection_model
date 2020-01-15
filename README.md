# Coinfection model

Illustrates how to share data between two unrelated models.

## Structure

* `hiv_model.py` contains the "HIV" model (actually not, since it's SIR).
* `tb_model.py` contains the "TB" model (actually not, since it's SI).
* `shared_vars.py` defines the variables that will be shared between the models.
* `run.py` coordinates the running of the models either independently or together.

## Usage

`python -i run.py`