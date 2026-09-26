/**
 * The radial tree of life: ~560 taxa laid out by `d3-hierarchy`.
 *
 * This is a controlled, presentational component. It owns no filter box, no
 * selection state and no navigation — `TreeExplorerPage` owns all three and
 * every control is real HTML on that page, never an SVG overlay that fades in
 * on hover. The only outward calls are `onSelect` and `onMatchStats`.
 */

import { useCallback, useEffect, useImperativeHandle, useRef } from 'react';
import type { Ref } from 'react';
import { select } from 'd3-selection';
import { hierarchy, tree as d3tree } from 'd3-hierarchy';
import { zoom as d3zoom, zoomIdentity } from 'd3-zoom';
import { arc as d3arc, linkRadial } from 'd3-shape';
import { scaleOrdinal } from 'd3-scale';
import { hsl } from 'd3-color';
import 'd3-transition'; // side effect: Selection.prototype.transition
import type { HierarchyPointLink, HierarchyPointNode } from 'd3-hierarchy';
import type { D3ZoomEvent, ZoomBehavior } from 'd3-zoom';

import { treeData } from './taxa';
import { readTreeTheme } from './treeTheme';
import { cn } from '../../utils/cn';
import type { TaxonomicLevel, TreeNode } from '../../types/taxonomy';

type Node = HierarchyPointNode<TreeNode>;
type Link = HierarchyPointLink<TreeNode>;

/** Above this, 14px species names overplot into mush; render none instead. */
const MAX_SPECIES_LABELS = 80;
/** Ring clearance beyond the outermost tip, sized for the new type scale. */
const LABEL_MARGIN = 140;
/** Below this the layout is unreadable anyway; stop shrinking. */
const MIN_RADIUS = 80;
/** Clearance from the outermost tip to the ring of group labels. */
const GROUP_LABEL_CLEARANCE = 22;
/** A group needs this many species before it earns a ring label and an arc. */
const MIN_GROUP_SIZE = 10;
/** Species tick offset, one step wider than the 9px original's 12. */
const TIP_LABEL_OFFSET = 14;

const LINK_WIDTH_MATCH = 1.2;
const LINK_WIDTH_BASE = 0.6;
const LINK_OPACITY_MUTED = 0.2;

/** Group arc geometry. Typed so the arc generator needs no cast. */
interface GroupArc {
  groupName: string;
  startAngle: number;
  endAngle: number;
  /** Ring radius: the arc has zero thickness, so inner === outer === this. */
  radius: number;
  color: string;
  nodeCount: number;
  hasMatches: boolean;
}

interface ArcEndpoint extends GroupArc {
  angle: number;
}

interface GroupPosition {
  groupName: string;
  hasMatches: boolean;
  /** Percentage along the circular text path. */
  pathOffset: number;
}

export interface TreeHandle {
  zoomIn: () => void;
  zoomOut: () => void;
  reset: () => void;
}

export interface TreeProps {
  /** Case-insensitive substring filter over species names. '' shows all. */
  filter: string;
  /** Rank used for grouping, colouring, and the ring arcs. */
  level: TaxonomicLevel;
  /** 'inspect' enables hover/click on tips; 'pan' enables drag-to-pan. */
  interaction: 'inspect' | 'pan';
  selected: string | null;
  onSelect: (species: string | null) => void;
  onMatchStats: (count: number, suppressed: boolean) => void;
  ref?: Ref<TreeHandle>;
}

/** Evenly spaced hues, one per group, at a fixed saturation and lightness. */
function generateColors(count: number): string[] {
  return Array.from({ length: count }, (_, i) => {
    const hue = ((i * 360) / count) % 360;
    return `${hsl(hue, 0.725, 0.475).formatHex()}bb`;
  });
}

export function Tree({
  filter,
  level,
  interaction,
  selected,
  onSelect,
  onMatchStats,
  ref: handleRef,
}: TreeProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const zoomRef = useRef<ZoomBehavior<SVGSVGElement, unknown> | null>(null);
  // The live view. A redraw (filter, level, selection, resize) rebuilds the SVG
  // from scratch, so without this the user's pan and zoom would be thrown away
  // every time they clicked a tip. Only the Reset control returns it to identity.
  const transformRef = useRef(zoomIdentity);

  // Held in a ref, not a dependency: the page passes an inline arrow, and a
  // draw that calls it would otherwise invalidate `drawTree` and redraw forever.
  const onMatchStatsRef = useRef(onMatchStats);
  useEffect(() => {
    onMatchStatsRef.current = onMatchStats;
  }, [onMatchStats]);

  const drawTree = useCallback(() => {
    if (!containerRef.current || !svgRef.current) return;

    const { width, height } = containerRef.current.getBoundingClientRect();
    // jsdom has no layout, and the first frame is measured before paint.
    if (width === 0 || height === 0) return;

    const theme = readTreeTheme(containerRef.current);
    const radius = Math.max(MIN_RADIUS, Math.min(width, height) / 2 - LABEL_MARGIN);

    const svg = select(svgRef.current);
    svg.selectAll('*').remove();

    // Centred viewBox: (0,0) is the middle of the canvas, so the layout group
    // needs no translate and an untouched view is exactly the identity transform.
    svg
      .attr('viewBox', `${-width / 2} ${-height / 2} ${width} ${height}`)
      .attr('width', width)
      .attr('height', height);

    const container = svg.append('g');

    const zoomBehavior = d3zoom<SVGSVGElement, unknown>()
      .scaleExtent([0.4, 12])
      .on('zoom', (event: D3ZoomEvent<SVGSVGElement, unknown>) => {
        transformRef.current = event.transform;
        container.attr('transform', event.transform.toString());
      })
      // The lock is gone. Pan/zoom is enabled exactly when the user asked for it.
      // Programmatic `zoom.transform` is unaffected, so the toolbar still works.
      .filter(() => interaction === 'pan');

    zoomRef.current = zoomBehavior;
    svg.call(zoomBehavior).call(zoomBehavior.transform, transformRef.current);

    const defs = svg.append('defs');
    const gradient = defs
      .append('radialGradient')
      .attr('id', 'rg-tree-background')
      .attr('cx', '50%')
      .attr('cy', '50%')
      .attr('r', '50%');
    gradient
      .append('stop')
      .attr('offset', '0%')
      .attr('stop-color', '#ffeaa7')
      .attr('stop-opacity', 0.33);
    gradient
      .append('stop')
      .attr('offset', '67%')
      .attr('stop-color', '#ffeaa7')
      .attr('stop-opacity', 0.11);
    gradient
      .append('stop')
      .attr('offset', '100%')
      .attr('stop-color', '#ffeaa7')
      .attr('stop-opacity', 0);

    container
      .append('circle')
      .attr('cx', 0)
      .attr('cy', 0)
      .attr('r', radius + 100)
      .style('fill', 'url(#rg-tree-background)')
      .style('pointer-events', 'none');

    const g = container.append('g');

    const layout = d3tree<TreeNode>()
      .size([2 * Math.PI, radius])
      .separation((a, b) => (a.parent === b.parent ? 1 : 2) / a.depth);
    const root = layout(hierarchy(treeData));

    /** The name of `node`'s ancestor (or self) at the grouping level. */
    const groupOf = (node: Node): string | null => {
      let current: Node | null = node;
      while (current) {
        if (current.data.taxonomicLevel === level) return current.data.name;
        current = current.parent;
      }
      return null;
    };

    const levelNames = new Set<string>();
    for (const node of root.descendants()) {
      if (node.data.taxonomicLevel === level) levelNames.add(node.data.name);
    }
    const colorScale = scaleOrdinal<string, string>()
      .domain([...levelNames])
      .range(generateColors(levelNames.size));

    // --- one match pass per draw, not one subtree walk per link -------------
    const q = filter.trim().toLowerCase();
    const matched = new Set<Node>();
    if (q !== '') {
      root.eachAfter((node: Node) => {
        const self =
          node.data.taxonomicLevel === 'species' && node.data.name.toLowerCase().includes(q);
        const kid = (node.children ?? []).some((child) => matched.has(child));
        if (self || kid) matched.add(node);
      });
    }
    const isMatch = (node: Node) => q === '' || matched.has(node);

    const speciesNodes = root.descendants().filter((d) => d.data.taxonomicLevel === 'species');
    const matching = q === '' ? [] : speciesNodes.filter((d) => matched.has(d));
    const suppressed = matching.length > MAX_SPECIES_LABELS;
    onMatchStatsRef.current(matching.length, suppressed);

    // --- links --------------------------------------------------------------
    const linkWidth = (link: Link) =>
      q !== '' && matched.has(link.target) ? LINK_WIDTH_MATCH : LINK_WIDTH_BASE;
    const linkOpacity = (link: Link) => (isMatch(link.target) ? 1 : LINK_OPACITY_MUTED);

    const linkPath = linkRadial<Link, Node>()
      .angle((d) => d.x)
      .radius((d) => d.y);

    const links = g
      .selectAll<SVGPathElement, Link>('.rg-tree__link')
      .data(root.links())
      .enter()
      .append('path')
      .attr('class', 'rg-tree__link')
      .attr('d', (d: Link) => linkPath(d) ?? '')
      .style('fill', 'none')
      .style('stroke', (d: Link) => {
        const group = groupOf(d.target);
        return group ? colorScale(group) : '#999999bb';
      })
      .style('stroke-width', linkWidth)
      .style('stroke-linecap', 'round')
      .style('opacity', linkOpacity);

    // --- species labels -----------------------------------------------------
    const tipTransform = (d: Node) => `rotate(${(d.x * 180) / Math.PI - 90}) translate(${d.y},0)`;
    const labelData = q === '' ? speciesNodes : suppressed ? [] : matching;

    const labels = g
      .selectAll<SVGTextElement, Node>('.rg-tree__species-label')
      .data(labelData)
      .enter()
      .append('text')
      .attr('class', 'rg-tree__species-label')
      .attr('dy', '0.31em')
      .attr('x', (d: Node) => (d.x < Math.PI === !d.children ? TIP_LABEL_OFFSET : -TIP_LABEL_OFFSET))
      .attr('text-anchor', (d: Node) => (d.x < Math.PI === !d.children ? 'start' : 'end'))
      .attr(
        'transform',
        (d: Node) =>
          `rotate(${(d.x * 180) / Math.PI - 90}) translate(${d.y},0)` +
          (d.x >= Math.PI ? ' rotate(180)' : ''),
      )
      .style('font-family', theme.fontFamily)
      .style('font-size', `${theme.speciesSize}px`)
      .style('font-weight', theme.speciesWeight)
      .style('fill', theme.text)
      .style('pointer-events', 'none')
      .style('opacity', q === '' ? 0 : 1)
      .text((d: Node) => d.data.name);

    // --- group rings --------------------------------------------------------
    const groupMembers = new Map<string, Node[]>();
    for (const node of speciesNodes) {
      const group = groupOf(node);
      if (!group) continue;
      const bucket = groupMembers.get(group);
      if (bucket) bucket.push(node);
      else groupMembers.set(group, [node]);
    }

    const labelRadius =
      Math.max(...root.descendants().map((node) => node.y)) + GROUP_LABEL_CLEARANCE;

    const textDefs = g.append('defs');
    textDefs
      .append('path')
      .attr('id', 'rg-tree-label-circle')
      .attr(
        'd',
        `M 0 ${-labelRadius} A ${labelRadius} ${labelRadius} 0 1 1 0 ${labelRadius} ` +
          `A ${labelRadius} ${labelRadius} 0 1 1 0 ${-labelRadius}`,
      )
      .style('fill', 'none')
      .style('stroke', 'none');

    const groupPositions: GroupPosition[] = [];
    groupMembers.forEach((nodes, groupName) => {
      if (nodes.length <= MIN_GROUP_SIZE) return;
      // Place the label at the group's mean angle, on its closest member.
      const avgAngle = nodes.reduce((sum, node) => sum + node.x, 0) / nodes.length;
      const representative = nodes.reduce((closest, node) =>
        Math.abs(node.x - avgAngle) < Math.abs(closest.x - avgAngle) ? node : closest,
      );
      let angle = representative.x;
      if (angle > 2 * Math.PI) angle -= 2 * Math.PI;
      if (angle < 0) angle += 2 * Math.PI;
      groupPositions.push({
        groupName,
        hasMatches: q === '' || nodes.some((node) => matched.has(node)),
        pathOffset: (angle / (2 * Math.PI)) * 100,
      });
    });

    g.selectAll<SVGTextElement, GroupPosition>('.rg-tree__group-label--muted')
      .data(groupPositions.filter((d) => q !== '' && !d.hasMatches))
      .enter()
      .append('text')
      .attr('class', 'rg-tree__group-label rg-tree__group-label--muted')
      .style('font-family', theme.fontFamily)
      .style('font-size', `${theme.groupMutedSize}px`)
      .style('font-weight', theme.groupMutedWeight)
      .style('fill', theme.textMuted)
      .style('pointer-events', 'none')
      .append('textPath')
      .attr('href', '#rg-tree-label-circle')
      .attr('startOffset', (d: GroupPosition) => `${d.pathOffset}%`)
      .style('text-anchor', 'middle')
      .text((d: GroupPosition) => d.groupName);

    const groupLabels = g
      .selectAll<SVGTextElement, GroupPosition>('.rg-tree__group-label--active')
      .data(groupPositions.filter((d) => d.hasMatches))
      .enter()
      .append('text')
      .attr('class', 'rg-tree__group-label rg-tree__group-label--active')
      .style('font-family', theme.fontFamily)
      .style('font-size', `${theme.groupSize}px`)
      .style('font-weight', theme.groupWeight)
      .style('fill', (d: GroupPosition) => colorScale(d.groupName))
      .style('pointer-events', 'none');
    groupLabels
      .append('textPath')
      .attr('href', '#rg-tree-label-circle')
      .attr('startOffset', (d: GroupPosition) => `${d.pathOffset}%`)
      .style('text-anchor', 'middle')
      .text((d: GroupPosition) => d.groupName);

    // --- group arcs ---------------------------------------------------------
    const groupArcs: GroupArc[] = [];
    groupMembers.forEach((nodes, groupName) => {
      if (nodes.length <= MIN_GROUP_SIZE) return;
      const angles = nodes.map((node) => node.x).sort((a, b) => a - b);
      groupArcs.push({
        groupName,
        startAngle: angles[0],
        endAngle: angles[angles.length - 1],
        radius: Math.max(...nodes.map((node) => node.y)) + 7,
        color: colorScale(groupName),
        nodeCount: nodes.length,
        hasMatches: q === '' || nodes.some((node) => matched.has(node)),
      });
    });

    const arcOpacity = (d: GroupArc) => (q === '' ? 0.8 : d.hasMatches ? 1 : 0.3);

    const arcPath = d3arc<GroupArc>()
      .innerRadius((d) => d.radius)
      .outerRadius((d) => d.radius)
      .startAngle((d) => d.startAngle)
      .endAngle((d) => d.endAngle);

    g.selectAll<SVGPathElement, GroupArc>('.rg-tree__arc')
      .data(groupArcs)
      .enter()
      .append('path')
      .attr('class', 'rg-tree__arc')
      .attr('d', (d: GroupArc) => arcPath(d) ?? '')
      .style('fill', 'none')
      .style('stroke', (d: GroupArc) => d.color)
      .style('stroke-width', 2)
      .style('stroke-linecap', 'round')
      .style('opacity', arcOpacity)
      .append('title')
      .text((d: GroupArc) => `${d.groupName} (${d.nodeCount} species)`);

    /** The little tick that caps each end of a group arc. */
    function endpointMarkers(
      className: string,
      angleOf: (d: GroupArc) => number,
      offsetY: number,
    ) {
      const data: ArcEndpoint[] = groupArcs.map((a) => ({ ...a, angle: angleOf(a) }));
      g.selectAll<SVGRectElement, ArcEndpoint>(`.${className}`)
        .data(data)
        .enter()
        .append('rect')
        .attr('class', className)
        .attr('x', (d: ArcEndpoint) => (d.radius + 2) * Math.cos(d.angle - Math.PI / 2) - 1)
        .attr('y', (d: ArcEndpoint) => (d.radius + 2) * Math.sin(d.angle - Math.PI / 2) - offsetY)
        .attr('width', 4)
        .attr('height', 3)
        .attr('transform', (d: ArcEndpoint) => {
          const x = (d.radius + 2) * Math.cos(d.angle - Math.PI / 2);
          const y = (d.radius + 2) * Math.sin(d.angle - Math.PI / 2);
          return `rotate(${(d.angle * 180) / Math.PI - 90}, ${x}, ${y})`;
        })
        .style('fill', (d: ArcEndpoint) => d.color)
        .style('opacity', arcOpacity);
    }
    endpointMarkers('rg-tree__arc-start', (d) => d.startAngle, 0);
    endpointMarkers('rg-tree__arc-end', (d) => d.endAngle, 3);

    // --- hover behaviour, shared by the hit targets and the labels ----------
    function highlightPath(node: Node) {
      const onPath = new Set<Node>();
      let current: Node | null = node;
      while (current) {
        onPath.add(current);
        current = current.parent;
      }
      links
        .style('stroke-width', (link: Link) =>
          onPath.has(link.target) ? LINK_WIDTH_MATCH : linkWidth(link),
        )
        .style('opacity', (link: Link) => (onPath.has(link.target) ? 1 : LINK_OPACITY_MUTED));
    }

    function restoreLinks() {
      links.style('stroke-width', linkWidth).style('opacity', linkOpacity);
    }

    function revealLabel(node: Node) {
      labels.filter((d: Node) => d === node).style('opacity', 1);
      groupLabels.style('opacity', 0);
    }

    function restoreLabels() {
      labels.style('opacity', q === '' ? 0 : 1);
      groupLabels.style('opacity', 1);
    }

    // --- species hit targets: usable without typing --------------------------
    g.selectAll<SVGCircleElement, Node>('.rg-tree__hit')
      .data(speciesNodes)
      .enter()
      .append('circle')
      .attr('class', 'rg-tree__hit')
      .attr('r', 5)
      .attr('transform', tipTransform)
      .style('fill', 'transparent')
      .style('cursor', interaction === 'inspect' ? 'pointer' : 'inherit')
      .style('pointer-events', interaction === 'inspect' ? 'all' : 'none')
      .on('mouseenter', (_event: MouseEvent, d: Node) => {
        revealLabel(d);
        highlightPath(d);
      })
      .on('mouseleave', () => {
        restoreLabels();
        restoreLinks();
      })
      .on('click', (event: MouseEvent, d: Node) => {
        event.stopPropagation();
        onSelect(d.data.name);
      })
      .append('title')
      .text((d: Node) => d.data.name);

    // --- the persistent selection marker ------------------------------------
    g.selectAll<SVGCircleElement, Node>('.rg-tree__selected')
      .data(speciesNodes.filter((d) => d.data.name === selected))
      .enter()
      .append('circle')
      .attr('class', 'rg-tree__selected')
      .attr('r', 4)
      .attr('transform', tipTransform)
      .style('fill', 'none')
      .style('stroke', theme.accent)
      .style('stroke-width', 2)
      .style('pointer-events', 'none');
  }, [filter, level, interaction, selected, onSelect]);

  // A ResizeObserver, not a window listener: toggling full screen changes the
  // panel box without ever firing a window resize.
  useEffect(() => {
    const element = containerRef.current;
    if (!element) return;
    let frame = 0;
    const observer = new ResizeObserver(() => {
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(drawTree);
    });
    observer.observe(element);
    drawTree();
    return () => {
      cancelAnimationFrame(frame);
      observer.disconnect();
    };
  }, [drawTree]);

  const scaleBy = useCallback((k: number) => {
    if (!svgRef.current || !zoomRef.current) return;
    select(svgRef.current).transition().duration(250).call(zoomRef.current.scaleBy, k);
  }, []);

  useImperativeHandle(
    handleRef,
    () => ({
      zoomIn: () => scaleBy(1.5),
      zoomOut: () => scaleBy(1 / 1.5),
      reset: () => {
        if (!svgRef.current || !zoomRef.current) return;
        select(svgRef.current)
          .transition()
          .duration(400)
          .call(zoomRef.current.transform, zoomIdentity);
      },
    }),
    [scaleBy],
  );

  return (
    <div
      ref={containerRef}
      className={cn('rg-tree__canvas', interaction === 'pan' && 'rg-tree__canvas--pan')}
    >
      <svg ref={svgRef} />
    </div>
  );
}
