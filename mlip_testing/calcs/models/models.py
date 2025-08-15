"""Define classes for all models."""

# ruff: noqa: D101, D102, F401

from __future__ import annotations

import dataclasses
from pathlib import Path

import mlipx
from mlipx.nodes.generic_ase import Device
import yaml

ALL_MODELS = {}

# https://github.com/ACEsuit/mace
ALL_MODELS["MACE-MP-0"] = mlipx.GenericASECalculator(
    module="mace.calculators",
    class_name="mace_mp",
    device="auto",
    kwargs={"model": "medium"},
    # MLIPX-hub model path, adjust as needed
)

# https://github.com/MDIL-SNU/SevenNet
ALL_MODELS["7net-0"] = mlipx.GenericASECalculator(
    module="sevenn.sevennet_calculator",
    class_name="SevenNetCalculator",
    device="auto",
    kwargs={"model": "7net-0"},
)
ALL_MODELS["7net-mf-ompa-mpa"] = mlipx.GenericASECalculator(
    module="sevenn.sevennet_calculator",
    class_name="SevenNetCalculator",
    device="auto",
    kwargs={"model": "7net-mf-ompa", "modal": "mpa"},
)


# https://github.com/orbital-materials/orb-models
@dataclasses.dataclass
class OrbCalc:  # numpydoc ignore=GL08
    name: str
    device: Device | None = None
    kwargs: dict = dataclasses.field(default_factory=dict)

    def get_calculator(self, **kwargs):  # numpydoc ignore=GL08
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator
        import torch._dynamo

        torch._dynamo.config.suppress_errors = True
        torch._dynamo.disable()
        import os

        os.environ["TORCH_DISABLE_MODULE_HIERARCHY_TRACKING"] = "1"

        method = getattr(pretrained, self.name)
        if self.device is None:
            orbff = method(**self.kwargs)
            calc = ORBCalculator(orbff, **self.kwargs)
        elif self.device == Device.AUTO:
            orbff = method(device=Device.resolve_auto(), **self.kwargs)
            calc = ORBCalculator(orbff, device=Device.resolve_auto(), **self.kwargs)
        else:
            orbff = method(device=self.device, **self.kwargs)
            calc = ORBCalculator(orbff, device=self.device, **self.kwargs)
        return calc

    @property
    def available(self) -> bool:  # numpydoc ignore=GL08
        try:
            from orb_models.forcefield import pretrained
            from orb_models.forcefield.calculator import ORBCalculator

            return True
        except ImportError:
            return False


ALL_MODELS["orb-v2"] = OrbCalc(name="orb_v2", device="auto")
ALL_MODELS["orb-v3"] = OrbCalc(name="orb_v3_conservative_inf_omat", device="auto")

# https://github.com/CederGroupHub/chgnet
ALL_MODELS["chgnet"] = mlipx.GenericASECalculator(
    module="chgnet.model",
    class_name="CHGNetCalculator",
)
# https://github.com/microsoft/mattersim
ALL_MODELS["mattersim"] = mlipx.GenericASECalculator(
    module="mattersim.forcefield",
    class_name="MatterSimCalculator",
    device="auto",
)
# https://www.faccts.de/orca/
ALL_MODELS["orca"] = mlipx.OrcaSinglePoint(
    orcasimpleinput="PBE def2-TZVP TightSCF EnGrad",
    orcablocks="%pal nprocs 8 end",
    orca_shell="",
)

# https://gracemaker.readthedocs.io/en/latest/gracemaker/foundation/
ALL_MODELS["GRACE-2L-OMAT"] = mlipx.GenericASECalculator(
    module="tensorpotential.calculator",
    class_name="TPCalculator",
    device=None,
    kwargs={
        "model": "../../models/GRACE-2L-OMAT",
    },
    # MLIPX-hub model path, adjust as needed
)


@dataclasses.dataclass
class FairChemCalc:  # numpydoc ignore=GL08
    model_name: str
    task_name: str
    device: Device | str = "cpu"
    overrides: dict = dataclasses.field(default_factory=dict)

    def get_calculator(self):  # numpydoc ignore=GL08
        from fairchem.core import FAIRChemCalculator, pretrained_mlip
        # torch.serialization.add_safe_globals([slice])

        predictor = pretrained_mlip.get_predict_unit(
            self.model_name, device=self.device, overrides=self.overrides
        )
        return FAIRChemCalculator(predictor, task_name=self.task_name)

    @property
    def available(self) -> bool:  # numpydoc ignore=GL08
        try:
            from fairchem.core import pretrained_mlip

            return self.model_name in pretrained_mlip._MODEL_CKPTS.checkpoints
        except Exception:
            return False


ALL_MODELS["fairchem-uma-sm"] = FairChemCalc(
    model_name="uma-sm", task_name="omat", device="cpu"
)


@dataclasses.dataclass
class SumCalcWrapper:  # numpydoc ignore=GL08
    model_names: list[str]
    d3_kwargs: dict | None = None

    def get_calculator(self):  # numpydoc ignore=GL08
        from ase import units
        from ase.calculators.mixing import SumCalculator
        import torch
        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator

        calculators = []
        for model in self.model_names:
            if isinstance(model, str):
                calc = MODELS[model]
                if hasattr(calc, "get_calculator"):
                    calc = calc.get_calculator()
            elif isinstance(model, dict):
                class_name = model["class_name"]
                module = model.get("module", None)
                device = model.get("device", "auto")
                kwargs = model.get("kwargs", {})

                if class_name == "FAIRChemCalculator":
                    from fairchem.core import FAIRChemCalculator, pretrained_mlip

                    predictor = pretrained_mlip.get_predict_unit(
                        kwargs["model_name"],
                        device=device,
                        overrides=kwargs.get("overrides", {}),
                    )
                    calc = FAIRChemCalculator(
                        predictor, task_name=kwargs.get("task_name", "omat")
                    )

                elif class_name == "OrbCalc":
                    name = kwargs.pop("name")
                    calc = OrbCalc(
                        name=name, device=device, kwargs=kwargs
                    ).get_calculator()

                else:
                    calc = mlipx.GenericASECalculator(
                        module=module,
                        class_name=class_name,
                        device=device,
                        kwargs=kwargs,
                    )
                    if hasattr(calc, "get_calculator"):
                        calc = calc.get_calculator()
            else:
                raise TypeError(
                    f"Unsupported model type in SumCalculator: {type(model)}"
                )

            calculators.append(calc)

        # calculators = []
        # for name in self.model_names:
        #     calc = MODELS[name]
        #     if hasattr(calc, "get_calculator"):
        #         calc = calc.get_calculator()
        #     calculators.append(calc)

        if self.d3_kwargs:
            d3_calc = TorchDFTD3Calculator(
                device=self.d3_kwargs.get("device", "cpu"),
                damping=self.d3_kwargs.get("damping", "bj"),
                xc=self.d3_kwargs.get("xc", "pbe"),
                dtype=getattr(torch, self.d3_kwargs.get("dtype", "float32")),
                cutoff=self.d3_kwargs.get("cutoff", 40.0 * units.Bohr),
            )
            calculators.append(d3_calc)

        return SumCalculator(calculators)


# OPTIONAL
# ========
# If you have custom property names you can use the UpdatedFramesCalc
# to set the energy, force and isolated_energies keys mlipx expects.

# REFERENCE = mlipx.UpdateFramesCalc(
#     results_mapping={"energy": "DFT_ENERGY", "forces": "DFT_FORCES"},
#     info_mapping={mlipx.abc.ASEKeys.isolated_energies.value: "isol_ene"},
# )

# ============================================================
# THE SELECTED MODELS!
# ONLY THESE MODELS WILL BE USED IN THE RECIPE
# ============================================================
MODELS = {}

# Load models from registry YAML: models.yaml
with open(Path(__file__).parent / "models.yaml") as f:
    _registry = yaml.safe_load(f)

for _name, _cfg in _registry.items():
    print(f"Loading model from models.yaml: {_name}")

    if _cfg["class_name"] == "FAIRChemCalculator":
        kwargs = _cfg.get("kwargs", {})
        MODELS[_name] = FairChemCalc(
            model_name=kwargs["model_name"],
            task_name=kwargs.get("task_name", "omat"),
            device=_cfg.get("device", "cpu"),
            overrides=kwargs.get("overrides", {}),
        )
    elif _cfg["class_name"] == "OrbCalc":
        kwargs = _cfg.get("kwargs", {})
        MODELS[_name] = OrbCalc(
            name=kwargs["name"],
            device=_cfg.get("device", "cpu"),
        )
    elif _cfg["class_name"] == "SumCalculator":
        kwargs = _cfg.get("kwargs", {})
        MODELS[_name] = SumCalcWrapper(
            model_names=_cfg["models"],
            d3_kwargs=_cfg.get("d3_kwargs", None),
        )
    else:
        MODELS[_name] = mlipx.GenericASECalculator(
            module=_cfg["module"],
            class_name=_cfg["class_name"],
            device=_cfg.get("device", "auto"),
            kwargs=_cfg.get("kwargs", {}),
        )
