"""stdio entry point for local MCP usage."""

import argparse
import sys


def main(argv=None):
    # Parsed before the import below, so --help works without the 'mcp' extra.
    argparse.ArgumentParser(
        prog="refgenie-mcp",
        description=(
            "Serve your local refgenie database to AI assistants over MCP (stdio), "
            "read-only. It reads the same REFGENIE_* variables as the refgenie CLI; "
            "REFGENIE_HOME_PATH (default ~/.refgenie) picks the folder."
        ),
        epilog=(
            "Add it to Claude Code with refgenie's virtual environment active, so the\n"
            "full path is saved:\n"
            '  claude mcp add refgenie -e REFGENIE_HOME_PATH=$HOME/.refgenie -- "$(which refgenie-mcp)"\n'
            "The public catalog is at https://api.refgenie.org/mcp."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    ).parse_args(argv)

    # Imported lazily: the `mcp` package ships only in the 'mcp' extras, but this
    # console script is installed by every refgenie install. Report the missing
    # extras the way the main CLI does, instead of dumping an import traceback.
    try:
        from refgenie import Refgenie
        from refgenie.mcp import tools
    except ImportError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    tools.set_refgenie(Refgenie())
    # stderr: stdout carries JSON-RPC.
    print(f"refgenie-mcp: {tools._server_line()}", file=sys.stderr)
    tools.mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
