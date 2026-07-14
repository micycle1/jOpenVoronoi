[![](https://jitpack.io/v/micycle1/jOpenVoronoi.svg)](https://jitpack.io/#micycle1/jOpenVoronoi)

# JOpenVoronoi+

A fork of *Rogach*'s [port](https://github.com/Rogach/jopenvoronoi) of the original C++ library [openvoronoi](https://github.com/aewallin/openvoronoi).

- Incremental Voronoi diagrams supporting point and line-segment sites
- Linear, parabolic, and conic Voronoi bisectors, with adaptive curve sampling and arc-length calculation
- Voronoi diagram generators: labyrinths and L-systems
- Nearest-face and nearest-*N* face queries via a KD-tree
- Diagram filtering:
  - polygon-interior filtering for oriented polygon boundaries
  - configurable medial-axis branch pruning
- Medial-axis extraction for polygons and line-site arrangements
- Voronoi-based offset generation, producing closed offset loops composed of line and arc elements
- Pocket-machining utilities (`jopenvoronoi-pocket` module):
  - medial-axis pocket decomposition and traversal
  - maximum-inscribed-circle (MIC) sampling along medial-axis edges
  - continuous spiral toolpath generation for clearing polygonal pockets
- JTS integration (`jopenvoronoi-jts` module):
  - build a Voronoi diagram directly from a JTS `Polygon` or `MultiPolygon`
  - export Voronoi cells and interior cells as JTS geometries
  - extract, dissolve, and optionally prune medial-axis branches
  - generate medial-axis polygon coverage
- Convenience aggregate artifact (`jopenvoronoi-all`) that brings in the core, generator, pocket, and JTS modules

# Fork Changes

*Rogach*'s work laid a strong foundation... This fork productionises the original codebase, transforming the proof-of-concept port into a modern, high-performance library. The project has undergone a significant overhaul to improve codebase health, performance, and ease of use.

# Example code


```java
import org.rogach.jopenvoronoi.*;

VoronoiDiagram voronoi = new VoronoiDiagram();
for (int i = 0; i < 100; i++) {
  voronoi.insertPointSite(Math.random(), Math.random());
}

voronoi.getFaces().forEach(face -> {
    face.getEdges().forEach(edge -> {
        vertex(edge.source.position.x, edge.source.position.y);
        vertex(edge.target.position.x, edge.target.position.y);
    });
});
```

Get the medial axis for a polygon represented by a list of points:

```java
List<Point> polygon = List.of(
    new Point(-1.0, -1.0),
    new Point(1.0, -1.0),
    new Point(1.0, 1.0),
    new Point(-1.0, 1.0)
);

VoronoiDiagram voronoi = new VoronoiDiagram();
List<Vertex> polygonVertices = new ArrayList<>(polygon.size());
for (Point point : polygon) {
    polygonVertices.add(voronoi.insertPointSite(point));
}
for (int i = 0; i < polygonVertices.size(); i++) {
    voronoi.insertLineSite(polygonVertices.get(i), polygonVertices.get((i + 1) % polygonVertices.size()));
}

// The example polygon points are ordered counter-clockwise.
voronoi.filter(new PolygonInteriorFilter(true));
voronoi.filter(new MedialAxisFilter());

List<Edge> medialAxis = voronoi.getDiagram().edges.stream()
    .filter(edge -> edge.valid
        && edge.type != EdgeType.LINESITE
        && edge.type != EdgeType.NULLEDGE
        && edge.type != EdgeType.OUTEDGE)
    .collect(Collectors.toList());
```

`medialAxis` contains the half-edges of the polygon's medial axis. Each geometric segment appears twice in the half-edge diagram (once per direction), so pair `edge`/`edge.twin` if you only need one copy of each branch.

Get an inward offset polygon from the same point-list input:

```java
List<Point> polygon = List.of(
    new Point(-1.0, -1.0),
    new Point(1.0, -1.0),
    new Point(1.0, 1.0),
    new Point(-1.0, 1.0)
);

VoronoiDiagram voronoi = new VoronoiDiagram();
List<Vertex> polygonVertices = new ArrayList<>(polygon.size());
for (Point point : polygon) {
    polygonVertices.add(voronoi.insertPointSite(point));
}
for (int i = 0; i < polygonVertices.size(); i++) {
    voronoi.insertLineSite(polygonVertices.get(i), polygonVertices.get((i + 1) % polygonVertices.size()));
}

// The example polygon points are ordered counter-clockwise.
voronoi.filter(new PolygonInteriorFilter(true));

List<OffsetLoop> insetLoops = new Offset(voronoi.getDiagram()).offset(0.2);
```

`insetLoops` contains closed offset loops. Omit the `PolygonInteriorFilter` step if you want both the inward and outward offsets from the full diagram.

## JTS integration

The `jopenvoronoi-jts` module bridges jOpenVoronoi with [JTS](https://github.com/locationtech/jts). `JtsVoronoi` wraps a `Polygon`/`MultiPolygon` and exposes the common derived products directly as JTS geometries, handling site insertion, boundary-side filtering, and curved-edge sampling for you:

```java
import org.locationtech.jts.geom.*;
import org.rogach.jopenvoronoi.jts.JtsVoronoi;

Polygon polygon = ...; // any JTS Polygon or MultiPolygon

JtsVoronoi jv = new JtsVoronoi(polygon);

List<Polygon> cells         = jv.getCells();              // every Voronoi cell
List<Polygon> interiorCells = jv.getInteriorCells();       // cells generated by the polygon interior
MultiLineString axis        = jv.getMedialAxis();          // full medial axis, dissolved into branches
MultiLineString prunedAxis  = jv.getMedialAxis(0.5);       // low-salience branches pruned (0 = none, 1 = max)
List<Polygon> coverage      = jv.getMedialAxisCoverage();  // polygon cut along its medial axis
```

`getVoronoiDiagram()` returns the underlying `VoronoiDiagram` for lower-level access, and `JtsGeometryIO` offers the lower-level conversion primitives (site insertion from arbitrary geometries, `MultiLineString` export) that `JtsVoronoi` is built on.

# Images

## Voronoi

<table>
  <tr>
    <td align="center" valign="center">Voronoi from poisson disc points</td>
     <td align="center" valign="center">Voronoi from poisson disc points (bounded)</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/voronoi2.png"></td>
    <td valign="top"><img src="resources/voronoi1.png"></td>
  </tr>
 </table>

 <table>
  <tr>
    <td align="center" valign="center">Voronoi from random points</td>
     <td align="center" valign="center">3 line sites</td>
  </tr>
  <tr>
    <td valign="top" width="50%"><img src="resources/voronoi3.png"></td>
    <td valign="top" width="50%"><img src="resources/rotating.gif"></td>
  </tr>
 </table>

## Generators

<table>
  <tr>
    <td align="center" valign="center">Gosper Curve</td>
     <td align="center" valign="center">Moore Curve</td>
  </tr>
  <tr>
    <td valign="top"><img src="resources/lindenmayer.png" width=500></td>
    <td valign="top"><img src="resources/moore.png"></td>
  </tr>
 </table>

 <table>
  <tr>
    <td align="center" valign="center">Labryinth</td>
     <td align="center" valign="center">Medial Axis (green)</td>
  </tr>
  <tr>
    <td valign="top" width="50%"><img src="resources/labryinth.png"></td>
    <td valign="top" width="50%"><img src="resources/medialAxis.png"></td>
  </tr>
 </table>

<!-- <p float="middle">
  <img src="resources/lindenmayer.png" alt="" width="49%"/>
  <img src="resources/labryinth.png" alt="" width="49%"/>
  <img src="resources/moore.png" alt="" width="49%"/>
  <img src="resources/medialAxis.png" alt="" width="49%"/>
</p> -->

## Nearest Faces 

<p float="middle">
  <img src="resources/neighbours/single.gif" alt="" width="49%"/>
  <img src="resources/neighbours/multi.gif" alt="" width="49%"/>
</p>

### Known Limitations
* The libary does not support >2 line sites terminating at the same vertex. Such arrangements *can* work, but it is not guaranteed.

License
=======
JOpenVoronoi is released under GPLv3, just like its parent
 [openvoronoi](https://github.com/aewallin/openvoronoi) project.
