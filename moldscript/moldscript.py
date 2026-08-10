import subprocess, sys
from moldscript.files import files
from moldscript.fukui import fukui
from moldscript.opt import opt
from moldscript.spc import spc
from moldscript.nmr import nmr
from moldscript.nbo import nbo
from moldscript.substructure import substructure
from moldscript.get_df import get_df
from moldscript.min_max import min_max
from moldscript.sterics import sterics
from moldscript.ensemble import ensemble
from moldscript.charges import charges
from moldscript.lowe import lowe
from moldscript.MLIP import mlip
from moldscript.argument_parser import (
    command_line_args,
    moldscript_version,
    moldscript_ref,
    time_run,
    write_arguments_file,
)
from moldscript.boltz import boltz
from moldscript.fmo import fmo
import time
from moldscript.utils import (
    append_run_log,
    emit,
    initialize_run_log,
    molecule_keys,
    print_run_header,
    terminal_error,
    terminal_success,
    terminal_warning,
)

header = """
   • ▌ ▄ ·.       ▄▄▌  ·▄▄▄▄  .▄▄ ·  ▄▄· ▄▄▄  ▪   ▄▄▄·▄▄▄▄▄
   ·██ ▐███▪▪     ██•  ██▪ ██ ▐█ ▀. ▐█ ▌▪▀▄ █·██ ▐█ ▄█•██  
   ▐█ ▌▐▌▐█· ▄█▀▄ ██▪  ▐█· ▐█▌▄▀▀▀█▄██ ▄▄▐▀▀▄ ▐█· ██▀· ▐█.▪
   ██ ██▌▐█▌▐█▌.▐▌▐█▌▐▌██. ██ ▐█▄▪▐█▐███▌▐█•█▌▐█▌▐█▪·• ▐█▌·
   ▀▀  █▪▀▀▀ ▀█▄▀▪.▀▀▀ ▀▀▀▀▀•  ▀▀▀▀ ·▀▀▀ .▀  ▀▀▀▀.▀    ▀▀▀ 
                             Paton Research Group, CO 2024                              
"""

def checks():
    # this is a dummy import just to warn the user if Open babel is not installed
    try:
        command_run_1 = ["obabel", "-H"]
        subprocess.run(
            command_run_1, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        terminal_error(
            "x  Open Babel is not installed! You can install the program with 'conda install -c conda-forge openbabel'"
        )
        sys.exit()
    try:
        from rdkit.Chem import AllChem as Chem
    except ModuleNotFoundError:
        terminal_error(
            "x  RDKit is not installed! You can install the program with 'conda install -c conda-forge rdkit'"
        )



def main():
    # This chunk parses the CLI arguments and load user-defined arguments from command line
    args = command_line_args()
    tstart = time.time()

    data_dicts = {}

    initialize_run_log(args.output, moldscript_version, time_run, moldscript_ref, sys.argv[1:])
    print_run_header(moldscript_version, time_run, moldscript_ref, sys.argv[1:])
    arguments_file = write_arguments_file(args)
    if arguments_file:
        emit(f"Saved reproducible argument file to {arguments_file}", style="cyan")


    first_read = ""
    if args.link:
        # ALL DATA
        all_read = files(calc="link", path=args.link, program=args.program)
        if args.opt:
            opt_data = opt(all_read.file_data, program=args.program)
            data_dicts["opt"] = opt_data
        if args.nmr:
            nmr_data = nmr(all_read.file_data, program=args.program)
            data_dicts["nmr"] = nmr_data
        if args.nbo:
            nbo_data = nbo(all_read.file_data, program=args.program)
            data_dicts["nbo"] = nbo_data

    else:
        # OPT
        if args.opt:
            opt_read = files("opt", args.opt, data_dicts, args.suffix_opt)
            if first_read == '':
                first_read = opt_read.file_data
            opt_data = opt(opt_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = opt_data.file_data
        
        # SPC
        if args.spc:
            spc_read = files(calc="spc", path=args.spc, data_dict=data_dicts, suffix=args.suffix_spc)
            if first_read == '':
                first_read = spc_read.file_data
            spc_data = spc(spc_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = spc_data.file_data
        
        # Charges
        if args.charges:
            chg_read = files(calc="charges", path=args.charges, data_dict=data_dicts, suffix=args.suffix_charges)
            if first_read == '':
                first_read = chg_read.file_data
            chg_data = charges(chg_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = chg_data.file_data
        
        # FMO
        if args.fmo:
            fmo_read = files(calc="fmo", path=args.fmo, data_dict=data_dicts, suffix=args.suffix_fmo)
            if first_read == '':
                first_read = fmo_read.file_data
            fmo_data = fmo(fmo_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = fmo_data.file_data
        
        # NMR
        if args.nmr:
            nmr_read = files("nmr", args.nmr, data_dicts, args.suffix_nmr)
            if first_read == '':
                first_read = nmr_read.file_data
            nmr_data = nmr(nmr_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = nmr_data.file_data

        # NBO
        if args.nbo:
            nbo_read = files("nbo", args.nbo, data_dicts, args.suffix_nbo)
            if first_read == '':
                first_read = nbo_read.file_data
            nbo_data = nbo(nbo_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = nbo_data.file_data

        # MLIP / MACE-Polar extxyz
        if args.mlip_neutral:
            mlip_data = mlip(
                neutral=args.mlip_neutral,
                reduced=args.mlip_reduced,
                oxidized=args.mlip_oxidized,
                data_dict=data_dicts,
                output=args.output,
            )
            data_dicts = mlip_data.file_data

        # FUKUI
        if args.fukui_neutral and args.fukui_reduced and args.fukui_oxidized:
            emit(f"FUKUI paths: {[args.fukui_neutral, args.fukui_reduced, args.fukui_oxidized]}", style="cyan")
            fukui_read = files(calc="fukui", data_dict=data_dicts, path=[args.fukui_neutral, args.fukui_reduced, args.fukui_oxidized], suffix= [args.suffix_fukui_neutral, args.suffix_fukui_reduced, args.suffix_fukui_oxidized])
            fukui_data = fukui(fukui_read.file_data, data_dicts, output=args.output, workers=args.workers)
            data_dicts = fukui_data.data_dict

    ensemble_only = bool(args.ensemble) and not molecule_keys(data_dicts)
    if args.ensemble:
        try:
            ensemble_data = ensemble(
                path=args.ensemble,
                data_dict=data_dicts,
                radii=args.ensemble_radii,
                grid=args.ensemble_grid,
                temp=args.temp,
                include_h=args.ensemble_include_h,
                exclude=args.ensemble_exclude,
                suffix=args.suffix_ensemble,
                output=args.output,
                workers=args.workers,
            )
            data_dicts = ensemble_data.file_data
        except (OSError, ValueError) as exc:
            terminal_error(f"x  Could not analyze XYZ ensemble: {exc}")
            raise SystemExit(1) from exc

    if args.substructure != "" and molecule_keys(data_dicts) and first_read:
        substructure_read = files(data_dict=data_dicts, calc="substructure", path=args.opt, suffix=args.suffix_opt)
        data_dicts = substructure(substructure_read.file_data, data_dicts, args.substructure, output=args.output).file_data
    elif args.substructure != "" and molecule_keys(data_dicts):
        terminal_warning(
            "SMARTS substructure matching requires a quantum-chemistry "
            "structure; it is not applied to ensemble-only XYZ input."
        )

    if (
        (args.volume != False or args.vall != False)
        and molecule_keys(data_dicts)
        and first_read
    ):
        data_dicts = sterics(first_read, data_dicts, args.volume, args.vall, args.radius, output=args.output).dd
    elif (
        (args.volume != False or args.vall != False)
        and molecule_keys(data_dicts)
    ):
        terminal_warning(
            "--volume and --vall require quantum-chemistry structure files; "
            "use --ensemble_radii for ensemble XYZ buried volumes."
        )

    df_getter = None
    if molecule_keys(data_dicts):
        df_getter = get_df(
            data_dicts,
            substructure=args.substructure,
            prefix=args.output,
            bond_filter=args.no_bond_filter,
            no_mol=args.no_mol,
            no_atom=args.no_atom,
            no_bond=args.no_bond or ensemble_only,
            mol_vector=args.mol_vector,
        )
    else:
        terminal_error(
            "x  No supported calculation or ensemble input files were requested"
        )
        raise SystemExit(1)

    if ensemble_only and (args.boltz or args.min_max or args.lowe):
        terminal_warning(
            "--boltz, --min_max, and --lowe operate on quantum-chemistry "
            "conformer tables and are ignored for ensemble-only XYZ input."
        )
    else:
        if args.boltz and df_getter is not None:
            boltz(
                temp=args.temp,
                prefix=args.output,
                energies=df_getter.energies,
            )
        if args.min_max and df_getter is not None:
            min_max(
                temp=args.temp,
                cut=args.cut,
                prefix=args.output,
                energies=df_getter.energies,
            )
        if args.lowe and df_getter is not None:
            lowe(prefix=args.output, energies=df_getter.energies)
    
    tfin = time.time()
    message = f"MolDscript finished running in {round(tfin - tstart, 2)} seconds"
    append_run_log(message)
    terminal_success(message)
if __name__ == "__main__":
    checks()
    main()
