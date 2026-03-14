[![](https://jitpack.io/v/micycle1/jOpenVoronoi.svg)](https://jitpack.io/#micycle1/jOpenVoronoi)

# JOpenVoronoi+

A fork of *Rogach*'s [port](https://github.com/Rogach/jopenvoronoi) of the original C++ library [openvoronoi](https://github.com/aewallin/openvoronoi).

# Features

- Incremental Voronoi Point & Line Site Diagrams
- Voronoi Diagram Generators: Labryinth, L-Systems
- Diagram Medial Axis
- Nearest N cells to point

# Fork Changes

This fork has the following changes:

<details>
  <summary>Click to expand!</summary>
  
- Converts source code to Java 11
- Introduces cell-point methods:
  - `nearestFace()` 
  - `nearestFaces()`
- Splits the library into appropriate sub-packages
- Converts source comments into proper Javadoc comments (only the most important comments converted so far...)
- Removes the constraint that point sites had to be placed within a unit-circle centered on (0,0) — now points can have any coordinate! (this could have side-effects...)
- Introduces *LindenmayerCurve*, *RandomLabyrinth* and *RandomPolygon* diagram generators (from the original's tests) into the main library under the *generate* sub-package
- Introduces `buildIntoVoronoiDiagram()` for PlanarGraphs
- Adds Javadoc comments to important arguments on generator classes
- Removes SVG output functionality
- Removes the debugging `step` argument (that was left in the code) from the main point/site insert methods
- Implements `position()` on `Edge`, `LineSite` and `Pointsite` classes
- Replace diagram's `HashSets` with `ArrayLists` for easier iteration
- More error handling
</details>



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

License
=======
JOpenVoronoi is released under GPLv3, just like its parent
 [openvoronoi](https://github.com/aewallin/openvoronoi) project.
