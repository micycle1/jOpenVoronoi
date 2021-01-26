# JOpenVoronoi+

A fork of *Rogach*'s [port](https://github.com/Rogach/jopenvoronoi) of the original C++ library [openvoronoi](https://github.com/aewallin/openvoronoi) — a library capable of creating 2D Voronoi segment diagrams and medial axes.

# Fork Changes

This fork has the following changes:

- Converts source code to Java 11
- Removes the *tests* maven sub-module and brings the main library up to the top level (so it's easily hostable as an artifact via JitPack)
- Splits the library into appropriate sub-packages
- Converts source comments into proper Javadoc comments (only the most important comments converted so far...)
- Removes the constraint that point sites had to be placed within a unit-circle centered on (0,0) — now points can have any coordinate!
- Introduces *LindenmayerCurve*, *RandomLabyrinth* and *RandomPolygon* diagram generators (from the original's tests) into the main library under the *generate* sub-package
- Adds Javadoc comments to important arguments on generator classes
- Removes SVG output functionality
- Removes the debugging `step` argument (that was left in the code) from the main point/site insert methods
- Implements `position()` on `LineSite` and `Pointsite` classes

Example code
============

```java
import org.rogach.jopenvoronoi.*;

VoronoiDiagram vd = new VoronoiDiagram();
Vertex v1 = vd.insert_point_site(new Point(-0.4,-0.2));
Vertex v2 = vd.insert_point_site(new Point(0,0.4));
Vertex v3 = vd.insert_point_site(new Point(0.4,-0.2));
vd.insert_line_site(v1, v2);
vd.insert_line_site(v2, v3);
vd.insert_line_site(v3, v1);
SvgOutput.output(vd, "test.svg");
```

License
=======
JOpenVoronoi is released under GPLv3 (see COPYING), same as it's parent
 [openvoronoi](https://github.com/aewallin/openvoronoi) project.