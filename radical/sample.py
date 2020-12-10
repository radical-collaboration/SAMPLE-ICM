#!/usr/bin/env python

from radical.entk import Pipeline, Stage, Task, AppManager
from scipy.stats import randint
from numpy.random import choice
from time import time
import os
import json
import argparse
import itertools
import pandas as pd
import pyDOE as doe


class Sample:
    def __init__(self, resource_reqs, model, analysis, param_space, init_params, n_samples, n_iterations, convergence, output):
        # Get/Set radical configuration attributes
        if os.environ.get("RADICAL_ENTK_VERBOSE") == None:
            os.environ["RADICAL_ENTK_REPORT"] = "True"

        self.hostname = os.environ.get("RMQ_HOSTNAME", "localhost")
        self.port = os.environ.get("RMQ_PORT", 5672)
        self.username = os.environ.get("RMQ_USERNAME")
        self.password = os.environ.get("RMQ_PASSWORD")

        self.param_space = os.path.join(os.getcwd(), 'params.csv')
        
        with open(param_space, 'r') as f:
            param_list = json.load(f)

        # Extract params from json file
        n_params = len(param_list)

        params = {}

        for i in range(n_params):
            params[param_list[i]["name"]] = []
            for j in range(param_list[i]["start"], param_list[i]["stop"] + 1, param_list[i]["step"]):
                params[param_list[i]["name"]].append(j)
        
        # compute all elements in parameter space
        ps = itertools.product(*params.values)


        # write parameter space to file
        with open(self.param_space, 'w') as f:
            f.write(",".join(params.keys()))
            f.write(os.linesep)
            f.write("\n".join([",".join(x for x in ps)]))

        self.resource_dict = json.load(resource_reqs)
        self.model = model
        self.analysis = analysis
        self.param_space = param_space
        self.n_samples = n_samples
        self.n_iterations = n_iterations
        self.convergence = convergence
        self.output = output
        self.select_file = init_params

        self.cur_dir = os.path.dirname(os.path.abspath(__file__))

    # has stop criteria been met?
    # naive stop condition for now
    def evaluate(self, cond, value):
        if cond == "convergence":
            with open("%s/selection.csv" % self.cur_dir, "r") as f:
                if len(f.readlines()) > 0 :
                    return True
                else:
                    return False

        elif cond == "iteration":
            return value

    # generate new set of tasks
    def generate(self, fldr_name):

        tasks = []
        for i in range(self.n_samples):
            t = Task()
            t.executable = self.model
            t.arguments = [fldr_name, i, self.select_file]

            tasks.append(t)

        return tasks

    # Adds selection task to pipeline
    def selection(self, ps_file, select_file):

        tasks = []

        t = Task()
        t.pre_exec = ['/bin/cp {0} {1}'.format(self.param_space, ps_file)]
        t.executable = self.analysis
        t.arguments = [ps_file, self.n_samples, select_file]
        t.download_output_data = [select_file]

        tasks.append(t)

        return tasks


    def run(self):

        # Get experiment time to save outputs to specific folder
        # To remove eventually and find a better solution
        fldr_name = os.path.join(self.output, str(int(time())))
        ps_file = os.path.join(fldr_name, os.path.basename(self.param_space))
        select_name = "selection.csv"
        select_file = os.path.join(fldr_name, select_name)

        i = 0
        # check if stop criteria has been met
        # otherwise run model with updated set of parameters

        value = i < self.n_iterations if self.n_iterations is not None
        cond = "convergence" if self.convergence is not None else "iterations"
        while self.evaluate(cond, value):

            # Setting up the pipeline
            # Create a Pipeline object
            p = Pipeline()

            # Create a Stage object
            s1 = Stage()

            # generate tasks
            tasks = self.generate(fldr_name)

            # add tasks to stage
            s1.add_tasks(tasks)

            # create selection stage
            s2 = Stage()

            # create selection task(s)
            selec_tasks = self.selection(ps_file, select_file)
            s2.add_tasks(selec_tasks)

            # Add Stage to the Pipeline
            p.add_stages([s1, s2])

            # Create Application Manager
            appman = AppManager(
                hostname=self.hostname,
                port=self.port,
                username=self.username,
                password=self.password,
            )

            # Assign the workflow as a set or list of Pipelines to the Application Manager
            appman.workflow = set([p])

            # Assign resource request description to the Application Manager
            appman.resource_desc = self.resource_dict

            # Run the Application Manager
            appman.run()
        
            # TODO find a better solution
            self.select_file = select_file 

            i += 1
            value = i < self.n_iterations if self.n_iterations is not None


def main():

    # Input argument parsing
    parser = argparse.ArgumentParser("Launch icemodel ensemble")
    parser.add_argument("model", type=str, help="location of the model executable")
    parser.add_argument("analysis", type=str, help="location of selection algorithm executable")
    parser.add_argument(
        "param_space",
        type=str,
        help="a json value containing the limits of the parameter space",
    )
    parser.add_argument(
        "init_params",
        type=str,
        help="csv file containing the initial parameter values"
    )
    parser.add_argument(
        "resource_reqs", type=argparse.FileType("r"), help="Resource requirements json",
    )
    parser.add_argument(
        "n_samples", type=int, help="number of samples to execute at a time."
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--n_iterations", type=int, help="maximum number of sampling iterations")
    group.add_argument("--convergence", action=store_true, help="execute until convergence")

    parser.add_argument(
        "output", type=str, help="output data location"
    )
    args = parser.parse_args()

    s = Sample(resource_reqs=args.resource_reqs, model=args.model, analysis=args.analysis, param_space=args.param_space,
            init_params=args.init_params, n_samples=args.n_samples, n_iterations=args.n_iterations, convergence=args.convergence, output=args.output)
    s.run()


if __name__ == "__main__":
    main()
