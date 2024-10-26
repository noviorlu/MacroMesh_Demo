
# README for CS 488/688 Assignment 3: Puppet

## Overview
This project implements a hierarchical, articulated 3D humanoid model, or "puppet," for CS 488/688 (Introduction to Computer Graphics) at the University of Waterloo. The puppet consists of various body parts connected by joints, allowing interactive control of the puppet's movements. The puppet is constructed using simple geometric primitives and transformed instances, arranged in a hierarchical data structure.

The application also includes interactive lighting, materials, and a Deferred Shading pipeline that supports multiple light sources, making the 3D structure of the puppet clear and visible without a significant performance hit.

### Features:
- **Hierarchical Model**: The puppet is built with a parent-child relationship between body parts.
- **Articulated Joints**: You can control at least 15 degrees of freedom (DOFs) across different joints.
- **Interactive Controls**: Manipulate the puppet with mouse and keyboard input.
- **Deferred Shading Pipeline**: Supports multiple light sources without a drop in framerate.
- **Basic Reset Function**: One "Reset All" function resets position, orientation, and joints together.

## How to Use

### 1. Running the Application
The executable requires a command-line argument specifying the Lua file for the puppet model. Use the following command in the terminal to run the application:

```bash
./clusterLOD <path_to_puppet_file>
```

For example, to run the sample puppet provided in `Assets/simpleScene.lua`, use:

```bash
./clusterLOD Assets/simpleScene.lua
```

### 2. Controls

#### Puppet Manipulation:
- **Left Mouse Button**: Drag to translate the puppet along the X and Y axes.
- **Middle Mouse Button**: Drag to move the puppet along the Z-axis (depth).
- **Right Mouse Button**: Rotate the puppet using a virtual trackball interface.

#### Joint Manipulation (J Mode):
- **Left Mouse Button**: Select or deselect joints.
- **Middle Mouse Button**: Rotate selected joints by dragging vertically (up/down).
- **Right Mouse Button**: Rotate the head when selected.

#### Keyboard Shortcuts:
- **P**: Switch to puppet manipulation mode (Position/Orientation).
- **J**: Switch to joint manipulation mode.
- **A**: Reset all (position, orientation, and joints).
- **C**: Toggle display of trackball circle.
- **Z**: Toggle Z-buffering.
- **B**: Toggle backface culling.
- **F**: Toggle frontface culling.
- **Q**: Quit the application.

### 3. Deferred Shading
The application uses a deferred shading pipeline that splits geometry and lighting into separate passes. This allows the use of multiple light sources without compromising performance. To adjust lighting, you can modify the Lua script for your puppet model and add more light sources or change their properties.

### 4. Known Limitations
- **Reset Functions**: The reset functionality only resets everything at once ("Reset All"). Separate resets for position, orientation, and joints are not implemented.
- **Undo/Redo**: Undo/redo stack for joint manipulations is not implemented.

## Puppet Structure
The puppet's hierarchical structure consists of a tree where:

- **Root Node**: Represents the torso and provides the transformation for the entire model.
- **Child Nodes**: Represent limbs, head, and other body parts, with each part connected to its parent. For example, the lower arm rotates around the elbow, which is a child node of the upper arm.

### Materials and Lighting
The puppet uses distinct materials for different body parts, which are applied using the `gr.material()` function in the Lua script. The Deferred Shading pipeline ensures that multiple light sources illuminate the puppet without causing significant frame drops.

## Additional Notes
- **Deferred Shading**: Implemented to enhance the visual quality of the puppet by handling multiple light sources.
- **Performance**: The Deferred Shading implementation avoids a significant framerate drop due to the efficient separation of geometry and lighting calculations.

## Conclusion
This application provides a basic framework for building and manipulating a 3D articulated puppet using Lua scripts. Although some features, such as undo/redo and separate reset functionalities, are not implemented, the application achieves good performance and supports dynamic lighting through a custom Deferred Shading pipeline.

For further extensions, you may add more complex joint structures or expand the lighting and material properties.

Enjoy exploring and manipulating your puppet!
