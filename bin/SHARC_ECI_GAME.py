#!/usr/bin/env python3
"""SHARC_ECI_GAME - overlap + phases extension to SHARC_ECI (Phase 1)."""

import os
import pickle
import datetime
import numpy as np

from SHARC_ECI import SHARC_ECI as _Base
import ECI as _ECI
from utils import phase_correction_cmplx

AUTHORS     = "eci-overlap-prototype"
VERSION     = "0.1-phase1"
VERSIONDATE = datetime.datetime(2026, 1, 1)
NAME        = "ECI_GAME"
DESCRIPTION = "SHARC_ECI subclass: overlap + phases (ECI2GAME Phase 1)"
CHANGELOG   = "Phase 1: overlap + phases extension. Basis-stable only."

_FEATURES      = {"h", "dm", "mol", "point_charges", "overlap", "phases"}
_PREV_CACHE    = {}
_CURRENT_SHARC = [None]


def _esd_spec(esd, sites):
    site_to_idx = {site: i for i, site in enumerate(sites)}
    spec = []
    for site, state in esd.site_states.items():
        si = site_to_idx[site]
        child = site.children_by_charge.get(state.Z)
        if child is None:
            raise RuntimeError(
                f"no child interface for fragment {site.label} charge {state.Z}"
            )
        idx = child.states.index(state)
        spec.append((si, idx, state.Z))
    return tuple(spec)


def _esd_overlap(prev_spec, curr_spec, sites):
    if len(prev_spec) != len(curr_spec):
        raise NotImplementedError("ESD fragment count changed")
    val = complex(1.0, 0.0)
    for (si_p, idx_p, z_p), (si_c, idx_c, z_c) in zip(prev_spec, curr_spec):
        if si_p != si_c or z_p != z_c:
            raise NotImplementedError("ESD site/charge ordering changed")
        child = sites[si_p].children_by_charge.get(z_p)
        if child is None or child.QMout["overlap"] is None:
            raise RuntimeError(
                f"child for fragment {sites[si_p].label} (charge {z_p}) "
                f"has no QMout['overlap']"
            )
        val *= complex(child.QMout["overlap"][idx_p, idx_c])
    return val


def _calculate_overlap(self):
    if self.previous is None:
        return
    for m in self.job.multiplicities:
        prev_specs = self.previous["esd_specs"][m]
        prev_psi   = self.previous["psi_esd"][m]
        curr_esds  = self.ECIbasis[m].ESDs
        if len(curr_esds) != len(prev_specs):
            raise NotImplementedError(
                f"ECI basis changed for mult {m}: "
                f"{len(prev_specs)} -> {len(curr_esds)}. "
                "Basis-changing overlap is Phase-2 work."
            )
        curr_specs = [_esd_spec(esd, self.sites) for esd in curr_esds]
        n = len(curr_esds)
        O_det = np.zeros((n, n), dtype=complex)
        for i, ps in enumerate(prev_specs):
            for j, cs in enumerate(curr_specs):
                O_det[i, j] = _esd_overlap(ps, cs, self.sites)
        curr_psi = self.ECIbasis[m].U @ self.Psi["ECSF"][m]
        self.overlap[m] = prev_psi.T @ O_det @ curr_psi


def _extract_state(job):
    state = {"esd_specs": {}, "psi_esd": {}}
    for m in job.job.multiplicities:
        esds = job.ECIbasis[m].ESDs
        state["esd_specs"][m] = [_esd_spec(esd, job.sites) for esd in esds]
        state["psi_esd"][m]   = job.ECIbasis[m].U @ job.Psi["ECSF"][m]
    return state


# ---------------------------------------------------------------- patches

_orig_eci_init = _ECI.ECI.__init__

def _patched_eci_init(self, *args, **kwargs):
    _orig_eci_init(self, *args, **kwargs)
    sharc = _CURRENT_SHARC[0]
    self.previous = _PREV_CACHE.get(self.job.charge)
    self.overlap  = {m: None for m in self.job.multiplicities}
    if sharc is not None:
        Z = self.job.charge
        ssc = { (k[0], k[1]): iface
                for k, iface in sharc._kindergarden["SSC"].items()
                if k[2] == Z }
        for site in self.sites:
            site.children_by_charge = {}
            for z in site.states.keys():
                site.children_by_charge[z] = ssc.get((site.label, z))

_ECI.ECI.__init__ = _patched_eci_init


_orig_eci_run = _ECI.ECI.run

def _patched_eci_run(self):
    self.configure()
    self.construct_excitonic_basis()
    self.allocate()
    self.calculate_Hmat()
    self.calculate_eigenstates()
    _calculate_overlap(self)
    self.calculate_properties()
    self.print_states()

_ECI.ECI.run = _patched_eci_run
_ECI.ECI.calculate_overlap_game = _calculate_overlap


# ---------------------------------------------------------------- interface

class SHARC_ECI_GAME(_Base):
    _version         = VERSION
    _versiondate     = VERSIONDATE
    _authors         = AUTHORS
    _name            = NAME
    _description     = DESCRIPTION
    _changelogstring = CHANGELOG

    @staticmethod
    def name():        return NAME
    @staticmethod
    def description(): return DESCRIPTION
    @staticmethod
    def about():       return f"{NAME}\n{DESCRIPTION}"

    def get_features(self, KEYSTROKES=None):
        return super().get_features(KEYSTROKES) | _FEATURES

    def read_requests(self, requests_file="QM.in"):
        super().read_requests(requests_file)
        if self.QMin.requests["phases"]:
            self.QMin.requests["overlap"] = True
        step = self.QMin.save["step"] or 0
        if step > 0:
            for child in self._kindergarden["SSC"].values():
                child.QMin.requests["overlap"] = True
        for request in ("dm",):
            for child in self._kindergarden["SSC"].values():
                child.QMin.requests[request] = self.QMin.requests[request]

    def run(self):
        savedir = self.QMin.save["savedir"]
        step    = self.QMin.save["step"] or 0
        _PREV_CACHE.clear()
        if step > 0:
            for Z in self.charges_to_do:
                path = os.path.join(savedir, f"eci_wf.{Z}.{step-1}.pkl")
                if os.path.isfile(path):
                    with open(path, "rb") as f:
                        _PREV_CACHE[Z] = pickle.load(f)
                    self.log.print(f"ECI_GAME: loaded {path}")
                else:
                    self.log.print(
                        f"ECI_GAME: warning - {path} missing; overlap will be null"
                    )
        _CURRENT_SHARC[0] = self
        try:
            super().run()
        finally:
            _CURRENT_SHARC[0] = None
        os.makedirs(savedir, exist_ok=True)
        for Z, job in self.ECIjobs.items():
            if job is None:
                continue
            path = os.path.join(savedir, f"eci_wf.{Z}.{step}.pkl")
            with open(path, "wb") as f:
                pickle.dump(_extract_state(job), f)
            self.log.print(f"ECI_GAME: wrote {path}")

    def getQMout(self):
        QMout = super().getQMout()
        if not (self.QMin.requests["overlap"] or self.QMin.requests["phases"]):
            return QMout
        for Z, job in self.ECIjobs.items():
            if job is None:
                continue
            for m, ovlp in job.overlap.items():
                if ovlp is None:
                    continue
                idxs = [i for i, s in enumerate(self.states)
                        if s.Z == Z and s.S == m - 1]
                for a, ia in enumerate(idxs):
                    for b, ib in enumerate(idxs):
                        QMout["overlap"][ia, ib] = ovlp[a, b]
        if self.QMin.requests["phases"]:
            QMout["phases"] = phase_correction_cmplx(
                QMout["overlap"], flag="simple"
            )
        return QMout


if __name__ == "__main__":
    SHARC_ECI_GAME().main()
