# Contributing to CoVigator pipeline

We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Develop with Github

We use github to host code, to track issues and feature requests, as well as accept pull requests.

## We Use [Github Flow](https://guides.github.com/introduction/flow/index.html), So All Code Changes Happen Through Pull Requests

Pull requests are the best way to propose changes to the codebase (we use [Github Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests.
3. If you've changed pipeline or the interface, update the documentation.
4. Ensure the test suite passes.
5. Issue that pull request!

## Repository overview

The repository is organized as follows:
- `main.nf` holds all the Nextflow code
- `environment-yml` defines all conda dependencies
- `nextflow.config` contains the nextflow configuration
- `bin` folder contains the Python code for the variant calling over the FASTA file. Any other custom code shall be added here.
- `reference` folder contains all reference files for SARS-CoV-2 which are used by default. To use CoVigator pipeline for a different organism these resources will need to be provided.

## Run the tests

Execute the following to run the tests: 
```
make
```

All test data is in the repository, but you will need Nextflow and conda installed.

These tests are automated in a private Gitlab instance, we may migrate this to GitHub CI system in the future.

## Any contributions you make will be under the MIT Software License

In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's [issues](https://github.com/briandk/transcriptase-atom/issues)

We use GitHub issues to track public bugs. Report a bug by [opening a new issue](); it's that easy!

## Write bug reports with detail, background, and sample code

[This is an example](http://stackoverflow.com/q/12488905/180626) of a bug report I wrote, and I think it's not a bad model. Here's [another example from Craig Hockenberry](http://www.openradar.me/11905408), an app developer whom I greatly respect.

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can. [My stackoverflow question](http://stackoverflow.com/q/12488905/180626) includes sample code that *anyone* with a base R setup can run to reproduce what I was seeing
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

## License

By contributing, you agree that your contributions will be licensed under its MIT License.

## References

This document was adapted from https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62
