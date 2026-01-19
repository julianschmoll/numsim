# ToDo Projekt

- [ ] Inflow und Outflow Boundary Conditions implementieren (Bene)
    - [ ] Inflow BC
    - [ ] Outflow BC
- [ ] Function als Geschwindigkeitsinput implementieren (für debugging purposes)

- [ ] Precice integration testen (koppeln mit solverdummy) (julian)
  - Dafür wrapper oder ein bisschen code architektur überlegen
- [ ] Herleitung der Kraftvektoren auf Festkörper (fabio)
    - Implementierung soll nur in eine Richtung kraft auswirken (orthogonal zur Wand, parallel vernachlässigen wir)
    - Was wenn Festkörper Rand berührt?

- [ ] Konfiguration überlegen von Tube
- [ ] Überlegung wie wir Zellen klassifizieren (fluid, solid, interface)
    - Fabian: evtl levelset methode anschauen
    - Julian: evtl anhand von Distanzfeldern arbeiten
    - Bene: evtl anhand von geometrischen überlegungen (bounding box etc)

- [ ] Calculix installieren
  - [ ] Calculix Config Files checken und 2d tube erstellen (dafür bedarf es noch keiner preCICE koppellung)
- [ ] Fenics installieren und testen

- [ ] Fluid-Solver und Structure Solver koppeln mit preCICE

- [ ] Plotting irgendwie (dicke von tube an stelle etc.), da mehr überlegen

- [ ] Präsentation vorbereiten
  - Schöne Bilder/Videos von Simulationen
  - Precice
  - Architektur Adapter
  - Herleitung Kraftvektoren
- 