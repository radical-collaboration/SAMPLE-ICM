#!/usr/bin/env python

from radical.entk import Pipeline, Stage, Task, AppManager
import os
import argparse
import json


def add_tasks(crevliq_range, clifvmax_range, stage, executable):

    for c in crevliq_range:
        for v in clifvmax_range:
            # Create a Task object which creates a file named 'output.txt' of size 1 MB
            t = Task()
            t.executable = executable
            t.arguments = [c, v]

            stage.add_tasks(t)


def main():

    if os.environ.get("RADICAL_ENTK_VERBOSE") == None:
        os.environ["RADICAL_ENTK_REPORT"] = "True"

    hostname = os.environ.get("RMQ_HOSTNAME", "localhost")
    port = os.environ.get("RMQ_PORT", 5672)
    username = os.environ.get("RMQ_USERNAME")
    password = os.environ.get("RMQ_PASSWORD")

    parser = argparse.ArgumentParser("Launch icemodel ensemble")
    parser.add_argument(
        "executable", type=str, help="location of the model executable"
    )
    parser.add_argument(
        "clifvmax_step", type=int, help="CLIFVMAX paramater range step value"
    )
    parser.add_argument(
        "crevliq_step", type=int, help="CREVLIQ paramater range step value"
    )
    parser.add_argument("param_space", type=int, help="The value of the 2d parameter space")
    parser.add_argument(
        "resource_reqs",
        type=argparse.FileType("r"),
        help="Resource requirements json",
    )
    args = parser.parse_args()

    clifvmax_range = range(0, args.clifvmax_step * args.param_space, args.clifvmax_step)
    crevliq_range = range(0, args.crevliq_step * args.param_space, args.crevliq_step)

    # Setting up the pipeline
    # Create a Pipeline object
    p = Pipeline()

    # Create a Stage object
    s1 = Stage()
    
    add_tasks(crevliq_range, clifvmax_range, s1, args.executable)

    # Add Stage to the Pipeline
    p.add_stages(s1)

    # Create Application Manager
    appman = AppManager(
        hostname=hostname, port=port, username=username, password=password
    )

    # Assign the workflow as a set or list of Pipelines to the Application Manager
    appman.workflow = set([p])

    # Create a dictionary to describe our resource request for XSEDE Stampede
    res_dict = json.load(args.resource_reqs)

    # Assign resource request description to the Application Manager
    appman.resource_desc = res_dict

    # Run the Application Manager
    appman.run()


if __name__ == "__main__":
    main()
