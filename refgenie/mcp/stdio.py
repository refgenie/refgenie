"""stdio entry point for local MCP usage."""

import sys


def main():
    # Imported lazily: the `mcp` package ships only in the 'mcp' extras, but this
    # console script is installed by every refgenie install. Report the missing
    # extras the way the main CLI does, instead of dumping an import traceback.
    try:
        from refgenie.mcp.tools import mcp
    except ImportError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
