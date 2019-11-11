#!/usr/bin/env python3
# -*- coding: utf-8 -*-

msg = (
"""
Works only in the Docker container!!!
- Start the container
    sudo docker start {Containername}
- Start an interactive terminal in the container
    sudo docker exec -it {Containername} bash
- Start the tests in the container terminal
    cd /
    pytest -vv --cov=pipeline /tests/
"""
)

# functions run in parallel do not count for line coverage
# test these functions here

# primerQC
#MFEprimer_template, MFEprimer_nontarget, MFEprimer_assembly, get_seq_fromDB
# make_templateDB, make_assemblyDB


if __name__ == "__main__":
    print(msg)