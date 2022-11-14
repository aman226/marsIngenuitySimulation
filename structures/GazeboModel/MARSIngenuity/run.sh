#!/bin/bash

echo "Welcome to M.A.R.S. Model.Please run the Simulink Model for Visualization!"
sleep 2
export SVGA_VGPU10=0
gazebo ingenuityMARS.world --verbose


