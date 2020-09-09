#!/usr/bin/env python

from radical.entk import Pipeline, Stage, Task, AppManager
from scipy.stats import randint
from numpy.random import choice
import os
import json
import argparse
import itertools
import pandas as pd
import pyDOE as doe


class Sample:
    def __init__(self, resource_reqs, model, param_space, n_samples, n_iterations):
        # Get/Set radical configuration attributes
        if os.environ.get("RADICAL_ENTK_VERBOSE") == None:
            os.environ["RADICAL_ENTK_REPORT"] = "True"

        self.hostname = os.environ.get("RMQ_HOSTNAME", "localhost")
        self.port = os.environ.get("RMQ_PORT", 5672)
        self.username = os.environ.get("RMQ_USERNAME")
        self.password = os.environ.get("RMQ_PASSWORD")

        self.resource_dict = json.load(resource_reqs)
        self.model = model
        self.param_space = param_space
        self.n_samples = n_samples
        self.n_iterations = n_iterations

        # first iteration
        self.lhs = None

    # has stop criteria been met?
    # naive stop condition for now
    def evaluate(self):
        self.n_iterations -= 1
        if self.n_iterations < 0:
            return True
        return False

    # which ensemble members require further analysis (detect high-gradient areas here)
    # update parameter space?
    def select(self, analysis=None):
        # random selection implementation for now
        if analysis is None:
            # return the rows that have not been previously executed
            if self.lhs is not None:
                self.param_space = self.param_space[~self.param_space.isin(self.lhs)]

            self.lhs = self.lhsampling()

    # generate new set of tasks
    def generate(self):

        tasks = []
        for i, row in lhs.iterrows():
            t = Task()
            t.executable = self.model
            t.arguments = [int(time()), *row.values.tolist()]

            tasks.append(t)

        return tasks

    def run(self):

        # check if stop criteria has been met
        # otherwise run model with updated set of parameters
        while not self.evaluate():

            # Setting up the pipeline
            # Create a Pipeline object
            p = Pipeline()

            # Create a Stage object
            s1 = Stage()

            # generate tasks
            tasks = self.generate()

            # add tasks to stage
            s1.add_tasks(tasks)

            # Add Stage to the Pipeline
            p.add_stages(s1)

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

    def lhsampling(self):

        if self.n_samples > len(self.param_space.index):
            self.n_samples = len(self.param_space.index)

        unit_lhs = pd.DataFrame(
            doe.lhs(
                n=len(self.param_space.columns),
                samples=self.n_samples,
                criterion="maximin",
            ),
            columns=self.param_space.columns,
        )
        lhs_idx = pd.DataFrame(index=unit_lhs.index)
        lhs = pd.DataFrame(index=unit_lhs.index)

        for i in unit_lhs.columns:
            lhs_idx[i] = (
                randint(0, len(self.param_space[i].unique()))
                .ppf(unit_lhs[i])
                .astype(int)
            )
            lhs[i] = lhs_idx[i].apply(lambda x: self.param_space[i].unique()[x])

        return lhs


def get_paramspace(params):
    p = json.load(params)
    ps = itertools.product(*list(p.values()))
    return pd.DataFrame(ps, columns=p.keys())


def main():

    # Input argument parsing
    parser = argparse.ArgumentParser("Launch icemodel ensemble")
    parser.add_argument("model", type=str, help="location of the model executable")
    parser.add_argument(
        "params",
        type=argparse.FileType("r"),
        help="json file containing all possible values of each parameter",
    )
    parser.add_argument(
        "resource_reqs", type=argparse.FileType("r"), help="Resource requirements json",
    )
    parser.add_argument(
        "n_samples", type=int, help="number of samples to execute at a time."
    )
    parser.add_argument(
        "n_iterations", type=int, help="maximum number of sampling iterations."
    )
    args = parser.parse_args()

    param_space = get_paramspace(args.params)

    s = Sample(resources_reqs=args.resource_reqs, model=args.model, param_space=self.param_space,
            n_samples=self.n_samples, n_iterations=self.n_iterations)
    s.run()


if __name__ == "__main__":
    main()
