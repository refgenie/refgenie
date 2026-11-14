// refgenie.org edge: redirects, then the static SPA and its files (llms.txt, robots.txt).
//
//   ui.refgenie.org/*, www.refgenie.org/*   -> https://refgenie.org/*       (301)
//   refgenie.org/{refgenie,refget,decisions,legacy}/*
//                                           -> https://docs.refgenie.org/*  (301)
//   everything else                         -> static assets (SPA)

const DOCS_SECTIONS = ["/refgenie", "/refget", "/decisions", "/legacy"];

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

    return env.ASSETS.fetch(request);
  },
};
