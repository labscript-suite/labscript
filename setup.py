from setuptools import setup

VERSION_SCHEME = {
    "version_scheme": "release-branch-semver",
    "local_scheme": "node-and-date",
}

setup(use_scm_version=VERSION_SCHEME)
