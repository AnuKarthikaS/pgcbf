"""Shared pytest fixtures: configuration and the certified systems are built once per session."""
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from pgcbf.certificates import build, certificate  # noqa: E402


def load(name):
    with open(os.path.join(ROOT, "config", f"{name}.json")) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def cfg():
    return load("default")


@pytest.fixture(scope="session", params=["default", "r1_passive"])
def system(request):
    c = load(request.param)
    s = build(c, train=True)
    return c, s, certificate(c, s)
