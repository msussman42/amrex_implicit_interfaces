#!/usr/bin/env python3

import argparse
import re
import shutil
import sys
from pathlib import Path


MAX_STEP_RE = re.compile(
    r"^\s*max_step\s*=\s*(\d+)\b"
)

KRYLOV_RE = re.compile(
    r"^\s*amr\.LSA_nsteps_krylov_subspace_method\s*=\s*(\d+)\b"
)


def read_parameters(input_file: Path):
    """
    Read:

        max_step = N
        amr.LSA_nsteps_krylov_subspace_method = k

    Commented lines are ignored.

    Returns:
        N, k
    """

    N = None
    k = None

    with input_file.open("r") as f:
        for line_number, line in enumerate(f, start=1):

            if line.lstrip().startswith("#"):
                continue

            match = MAX_STEP_RE.match(line)
            if match:
                N = int(match.group(1))

            match = KRYLOV_RE.match(line)
            if match:
                k = int(match.group(1))

    if N is None:
        raise ValueError(
            f"Could not find an active 'max_step = N' in:\n"
            f"    {input_file}"
        )

    if k is None:
        raise ValueError(
            "Could not find an active "
            "'amr.LSA_nsteps_krylov_subspace_method = k' in:\n"
            f"    {input_file}"
        )

    if k < 1:
        raise ValueError(
            "amr.LSA_nsteps_krylov_subspace_method must be >= 1"
        )

    return N, k


def insert_restart_line(
    filename: Path,
    N: int,
    i: int,
):
    """
    Insert a restart line into the generated file.

    If a commented-out amr.restart line exists, insert the new restart
    line immediately underneath it.

    Otherwise, append the new restart line to the end of the file.

    Example:

        #amr.restart = chk00000
        amr.restart=chk00500LSA00001
    """

    text = filename.read_text()
    lines = text.splitlines(keepends=True)

    restart_name = f"chk{N:05d}LSA{i:05d}"

    commented_restart_re = re.compile(
        r"^(?P<indent>\s*)#\s*amr\.restart\s*="
    )

    new_lines = []
    inserted = False

    for line in lines:

        new_lines.append(line)

        if inserted:
            continue

        line_without_newline = line.rstrip("\r\n")

        match = commented_restart_re.match(
            line_without_newline
        )

        if match:

            indent = match.group("indent")

            if line.endswith("\r\n"):
                newline = "\r\n"
            else:
                newline = "\n"

            new_lines.append(
                f"{indent}amr.restart={restart_name}{newline}"
            )

            inserted = True

    # If there is no commented-out amr.restart line,
    # append the new restart line at the end.
    if not inserted:

        if new_lines:
            if not (
                new_lines[-1].endswith("\n")
                or new_lines[-1].endswith("\r\n")
            ):
                new_lines[-1] += "\n"

        new_lines.append(
            f"\namr.restart={restart_name}\n"
        )

    filename.write_text("".join(new_lines))


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Create LSA input files from an AMReX input file."
        )
    )

    parser.add_argument(
        "--input_file",
        "--input-file",
        dest="input_file",
        default="inputs.growthrate.LSA",
        help=(
            "Path to the original input file. "
            "Default: ./inputs.growthrate.LSA"
        ),
    )

    parser.add_argument(
        "--output_folder",
        "--output-folder",
        dest="output_folder",
        default=None,
        help=(
            "Directory where generated files should be written. "
            "If omitted, generated files are placed in the same "
            "directory as the input file."
        ),
    )

    args = parser.parse_args()

    input_file = Path(
        args.input_file
    ).expanduser().resolve()

    if not input_file.exists():
        print(
            "\nERROR: Input file does not exist:\n"
            f"    {input_file}\n",
            file=sys.stderr,
        )
        return 1

    if not input_file.is_file():
        print(
            "\nERROR: --input_file must refer to a file:\n"
            f"    {input_file}\n",
            file=sys.stderr,
        )
        return 1

    # If no output folder is specified,
    # use the same directory as the input file.
    if args.output_folder is None:

        output_dir = input_file.parent

    else:

        output_dir = Path(
            args.output_folder
        ).expanduser().resolve()

        if output_dir.exists() and not output_dir.is_dir():
            print(
                "\nERROR: --output_folder must refer to a directory:\n"
                f"    {output_dir}\n",
                file=sys.stderr,
            )
            return 1

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

    print()
    print("------------------------------------------------------------")
    print("LSA input-file generator")
    print("------------------------------------------------------------")
    print()
    print(f"Input file       : {input_file}")
    print(f"Output directory : {output_dir}")
    print()

    try:

        N, k = read_parameters(input_file)

    except ValueError as exc:

        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )

        return 1

    print("Parameters found:")
    print()
    print(f"    max_step                              = {N}")
    print(
        f"    amr.LSA_nsteps_krylov_subspace_method = {k}"
    )
    print(f"    checkpoint prefix                     = chk")
    print()

    print(f"Will create {k} LSA input file(s).")
    print()

    base_name = input_file.name

    for i in range(1, k + 1):

        checkpoint_name = (
            f"chk{N:05d}LSA{i:05d}"
        )

        output_name = (
            f"{base_name}.{checkpoint_name}"
        )

        output_file = output_dir / output_name

        try:

            shutil.copy2(
                input_file,
                output_file,
            )

            insert_restart_line(
                output_file,
                N,
                i,
            )

        except Exception as exc:

            print(
                f"\nERROR while creating:\n"
                f"    {output_file}\n"
                f"{exc}\n",
                file=sys.stderr,
            )

            return 1

        print(
            f"[{i}/{k}] Created: {output_file}"
        )

    print()
    print("------------------------------------------------------------")
    print(f"Successfully created {k} file(s).")
    print("------------------------------------------------------------")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())