// refgenie.org edge: redirects and llms.txt, then the static SPA.
//
//   ui.refgenie.org/*, www.refgenie.org/*   -> https://refgenie.org/*       (301)
//   refgenie.org/{refgenie,refget,decisions,legacy}/*
//                                           -> https://docs.refgenie.org/*  (301)
//   refgenie.org/llms.txt                   -> index for AI agents
//   everything else                         -> static assets (SPA)

const DOCS_SECTIONS = ["/refgenie", "/refget", "/decisions", "/legacy"];

const LLMS_TXT = `# refgenie

> refgenie is a catalog of reference genomes and the assets built from them: FASTA files, aligner indexes and annotations. Find a genome by name or species, then download its assets or the command that pulls them. This site is the web UI; the API is at https://api.refgenie.org.

- [SKILL.md](https://refgenie.org/SKILL.md): how to find a genome and get its files with the HTTP API, step by step
- [OpenAPI schema](https://api.refgenie.org/openapi.json): every HTTP endpoint
- [MCP server](https://api.refgenie.org/mcp): Streamable HTTP, read-only tools; list tools return 20 rows a page

## Docs

- [Connect an AI assistant](https://docs.refgenie.org/refgenie/mcp/): setting up the MCP server
- [refgenie CLI skill](https://docs.refgenie.org/SKILL.md): driving the refgenie command-line tool
- [refgenie documentation](https://docs.refgenie.org/refgenie/)
`;

export default {
  async fetch(request, env) {
    const url = new URL(request.url);

    if (url.hostname === "ui.refgenie.org" || url.hostname === "www.refgenie.org") {
      url.hostname = "refgenie.org";
      return Response.redirect(url.toString(), 301);
    }

    const path = url.pathname;
    if (DOCS_SECTIONS.some((s) => path === s || path.startsWith(s + "/"))) {
      url.hostname = "docs.refgenie.org";
      return Response.redirect(url.toString(), 301);
    }

    if (path === "/llms.txt") {
      return new Response(LLMS_TXT, {
        headers: { "content-type": "text/plain; charset=utf-8" },
      });
    }

    return env.ASSETS.fetch(request);
  },
};
