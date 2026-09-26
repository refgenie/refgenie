/**
 * Project-level links: refgenie the project, not this instance of it.
 *
 * These are frontend constants ON PURPOSE, not the `links` block of
 * `/service-info`: the production server reports
 * `docs: https://refgenie.databio.org` — live, but the RETIRED readthedocs
 * site. A server operator has no business rewriting "where are the refgenie
 * docs": the answer is the same on every instance, so it lives here where it
 * can be checked, and `UiConfig` has no `links` field at all.
 *
 * Every URL here was fetched and returned 200. Two traps, both verified:
 *   - github.com/refgenie/refgenie1 is PRIVATE; the public repo is
 *     github.com/refgenie/refgenie (it receives refgenie1's code in the release
 *     swap). Never link refgenie1.
 *   - blob/master/CHANGELOG.md and blob/master/CONTRIBUTING.md 404 today,
 *     because that repo is still the legacy 0.x tree. Link repo roots and the
 *     docs site instead. LICENSE.txt is the one blob path that resolves both
 *     before and after the swap.
 */

export const PROJECT_LINKS = {
  docs: 'https://docs.refgenie.org',
  github: 'https://github.com/refgenie/refgenie',
  issues: 'https://github.com/refgenie/refgenie/issues',
  license: 'https://github.com/refgenie/refgenie/blob/master/LICENSE.txt',
  pypi: 'https://pypi.org/project/refgenie/',
  docsSource: 'https://github.com/refgenie/refgenie-docs',
  registry: 'https://github.com/refgenie/refgenie-registry',
  /** The instance-agnostic skill file: the CLI and Python API, not this API. */
  toolSkill: 'https://docs.refgenie.org/SKILL.md',
  useWithAi: 'https://docs.refgenie.org/refgenie/use_with_ai/',
} as const;

export interface LinkEntry {
  readonly href: string;
  readonly label: string;
  readonly blurb: string;
}

/** The doc map, in the order a newcomer needs it. */
export const DOC_GROUPS: ReadonlyArray<{
  readonly heading: string;
  readonly items: readonly LinkEntry[];
}> = [
  {
    heading: 'Start here',
    items: [
      {
        href: 'https://docs.refgenie.org/refgenie/',
        label: 'Introduction',
        blurb: 'What refgenie is and the problem it solves.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/install/',
        label: 'Install',
        blurb: 'pip, the optional extras, and what each one turns on.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/overview/',
        label: 'Overview',
        blurb: 'Genomes, asset classes, assets, recipes — how the pieces fit.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/whats_new/',
        label: "What's new in v1",
        blurb: 'What changed from the 0.x series.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/upgrade/',
        label: 'Upgrade from 0.x',
        blurb: 'Moving an existing refgenie config onto version 1.',
      },
    ],
  },
  {
    heading: 'Using refgenie',
    items: [
      {
        href: 'https://docs.refgenie.org/refgenie/pull/',
        label: 'Download pre-built assets',
        blurb: 'Pull from a public server onto your own machine.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/asset_registry_paths/',
        label: 'Refer to assets',
        blurb: 'Registry paths: the genome/asset_group.asset:tag syntax.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/seek/',
        label: 'Retrieve a path',
        blurb: 'Turn a registry path into a file path a pipeline can use.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/build/',
        label: 'Build assets',
        blurb: 'Run a recipe to produce an asset from its inputs.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/data_channels/',
        label: 'Use data channels',
        blurb: 'Where recipes and asset classes come from.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/mcp/',
        label: 'Connect AI assistants',
        blurb: 'The Model Context Protocol server.',
      },
    ],
  },
  {
    heading: 'Running an instance',
    items: [
      {
        href: 'https://docs.refgenie.org/refgenie/dash/',
        label: 'Use the dashboard',
        blurb: 'refgenie dash: this same interface over your own machine.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/server/',
        label: 'Run a server',
        blurb: 'refgenie serve: the public REST API and web UI.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/servers/',
        label: 'Public servers',
        blurb: 'Instances you can pull from without running anything.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/bridge/',
        label: 'Connect a local refgenie',
        blurb: 'How a public page talks to the dash on your own machine.',
      },
    ],
  },
  {
    heading: 'Reference',
    items: [
      {
        href: 'https://docs.refgenie.org/refgenie/glossary/',
        label: 'Glossary',
        blurb: 'Every term this interface uses, defined once.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/faq/',
        label: 'FAQ',
        blurb: 'The questions that come up most.',
      },
      {
        href: 'https://docs.refgenie.org/refgenie/contributing/',
        label: 'Contributing',
        blurb: 'How to propose a change.',
      },
      {
        href: 'https://docs.refgenie.org/legacy/',
        label: 'Legacy (pre-1.0) docs',
        blurb: 'Documentation for the retired 0.x series.',
      },
    ],
  },
];

export interface Citation {
  readonly authors: string;
  readonly year: string;
  readonly title: string;
  readonly venue: string;
  readonly doi: string;
  readonly note: string;
  /** One-line plain text for the copy button. */
  readonly plain: string;
}

/**
 * From refgenie-docs/docs/refgenie/manuscripts.md; author lists, titles, venues
 * and years re-confirmed against doi.org content negotiation.
 */
export const CITATIONS: readonly Citation[] = [
  {
    authors: 'Stolarczyk M, Reuter VP, Smith JP, Magee NE, Sheffield NC',
    year: '2020',
    title: 'Refgenie: a reference genome resource manager',
    venue: 'GigaScience 9(2)',
    doi: 'https://doi.org/10.1093/gigascience/giz149',
    note: 'The introductory publication. Cite this one if you cite only one.',
    plain:
      'Stolarczyk M, Reuter VP, Smith JP, Magee NE, Sheffield NC. ' +
      'Refgenie: a reference genome resource manager. GigaScience. 2020;9(2). ' +
      'doi:10.1093/gigascience/giz149',
  },
  {
    authors: 'Stolarczyk M, Xue B, Sheffield NC',
    year: '2021',
    title: 'Identity and compatibility of reference genome resources',
    venue: 'NAR Genomics and Bioinformatics 3(2)',
    doi: 'https://doi.org/10.1093/nargab/lqab036',
    note: 'Sequence-derived genome identifiers and provenance tracking — the work behind the digests this interface shows.',
    plain:
      'Stolarczyk M, Xue B, Sheffield NC. Identity and compatibility of ' +
      'reference genome resources. NAR Genomics and Bioinformatics. 2021;3(2). ' +
      'doi:10.1093/nargab/lqab036',
  },
];

export const ECOSYSTEM: readonly LinkEntry[] = [
  {
    href: 'https://docs.refgenie.org/refget/',
    label: 'refget',
    blurb:
      'The Python package behind refgenie genome identity: refget digests, sequence collections, and the RefgetStore format that holds the sequences themselves.',
  },
  {
    href: 'https://ga4gh.github.io/refget/sequences/',
    label: 'GA4GH refget sequences',
    blurb:
      'The standard for naming a sequence by the digest of its own content, so the same sequence gets the same identifier everywhere.',
  },
  {
    href: 'https://ga4gh.github.io/refget/seqcols/',
    label: 'GA4GH sequence collections',
    blurb:
      'The standard that gives a whole assembly a digest, and defines how to compare two assemblies for compatibility. Refgenie identifies every genome this way.',
  },
  {
    href: 'https://seqcolapi.databio.org/docs',
    label: 'Sequence Collections API',
    blurb:
      'A reference implementation of the sequence collections standard you can query directly.',
  },
  {
    href: 'https://github.com/refgenie/refgenie-registry',
    label: 'refgenie-registry',
    blurb:
      'The community data channel: the genome definitions, asset classes, and build recipes that refgenie syncs from.',
  },
];
