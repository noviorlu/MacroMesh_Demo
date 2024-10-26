-- humanoid_puppet_with_suzanne.lua
-- This script creates a humanoid skeleton with suzanne as the head and 17 joints using cube, sphere, and suzanne models.

root = gr.node('root')

-- Materials
red = gr.material({1.0, 0.0, 0.0}, {0.5, 0.4, 0.4}, 25)
blue = gr.material({0.5, 0.42, 0.0}, {0.4, 0.4, 0.5}, 25)
green = gr.material({0.5, 0.0, 0.5}, {0.4, 0.5, 0.4}, 25)
white = gr.material({0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, 25)

-- Torso
torso_joint = gr.joint('torso_joint', {-90, 0, 90}, {-90, 0, 90})
root:add_child(torso_joint)

torso = gr.mesh('cube', 'torso')
torso:scale(1.0, 2.0, 0.5)  -- Torso scale
torso:set_material(blue)
torso_joint:add_child(torso)

-- Shoulders
shoulder_joint = gr.joint('shoulder_joint', {-45, 0, 45}, {-45, 0, 45})
shoulder_joint:translate(0.0, 1.0, 0.0)
torso_joint:add_child(shoulder_joint)

shoulders = gr.mesh('cube', 'shoulders')
shoulders:scale(2.0, 0.5, 0.5)
shoulders:set_material(white)
shoulder_joint:add_child(shoulders)

-- Left Upper Arm
left_upper_arm_joint = gr.joint('left_upper_arm_joint', {-90, 0, 90}, {-90, 0, 90})
left_upper_arm_joint:translate(-1.0, 0.5, 0.0)
shoulder_joint:add_child(left_upper_arm_joint)

left_upper_arm = gr.mesh('cube', 'left_upper_arm')
left_upper_arm:scale(0.5, 1.5, 0.5)
left_upper_arm:translate(0.0, -1.5, 0.0)
left_upper_arm:set_material(red)
left_upper_arm_joint:add_child(left_upper_arm)

-- Left Forearm
left_forearm_joint = gr.joint('left_forearm_joint', {-90, 0, 90}, {-90, 0, 90})
left_forearm_joint:translate(0.0, -1.5, 0.0)
left_upper_arm_joint:add_child(left_forearm_joint)

left_forearm = gr.mesh('cube', 'left_forearm')
left_forearm:scale(0.5, 1.0, 0.5)
left_forearm:translate(0.0, -1.0, 0.0)
left_forearm:set_material(green)
left_forearm_joint:add_child(left_forearm)

-- Left Hand
left_hand_joint = gr.joint('left_hand_joint', {-45, 0, 45}, {-45, 0, 45})
left_hand_joint:translate(0.0, -1.0, 0.0)
left_forearm_joint:add_child(left_hand_joint)

left_hand = gr.mesh('cube', 'left_hand')
left_hand:scale(0.5, 0.5, 0.5)
left_hand:translate(0.0, -0.5, 0.0)
left_hand:set_material(blue)
left_hand_joint:add_child(left_hand)

-- Right Upper Arm
right_upper_arm_joint = gr.joint('right_upper_arm_joint', {-90, 0, 90}, {-90, 0, 90})
right_upper_arm_joint:translate(1.0, 0.5, 0.0)
shoulder_joint:add_child(right_upper_arm_joint)

right_upper_arm = gr.mesh('cube', 'right_upper_arm')
right_upper_arm:scale(0.5, 1.5, 0.5)
right_upper_arm:translate(0.0, -1.5, 0.0)
right_upper_arm:set_material(red)
right_upper_arm_joint:add_child(right_upper_arm)

-- Right Forearm
right_forearm_joint = gr.joint('right_forearm_joint', {-90, 0, 90}, {-90, 0, 90})
right_forearm_joint:translate(0.0, -1.5, 0.0)
right_upper_arm_joint:add_child(right_forearm_joint)

right_forearm = gr.mesh('cube', 'right_forearm')
right_forearm:scale(0.5, 1.0, 0.5)
right_forearm:translate(0.0, -1.0, 0.0)
right_forearm:set_material(green)
right_forearm_joint:add_child(right_forearm)

-- Right Hand
right_hand_joint = gr.joint('right_hand_joint', {-45, 0, 45}, {-45, 0, 45})
right_hand_joint:translate(0.0, -1.0, 0.0)
right_forearm_joint:add_child(right_hand_joint)

right_hand = gr.mesh('cube', 'right_hand')
right_hand:scale(0.5, 0.5, 0.5)
right_hand:translate(0.0, -0.5, 0.0)
right_hand:set_material(blue)
right_hand_joint:add_child(right_hand)

-- Hips
hips_joint = gr.joint('hips_joint', {-45, 0, 45}, {-45, 0, 45})
hips_joint:translate(0.0, -1.5, 0.0)
torso_joint:add_child(hips_joint)

hips = gr.mesh('cube', 'hips')
hips:scale(1.5, 0.5, 0.5)
hips:set_material(white)
hips_joint:add_child(hips)

-- Left Thigh
left_thigh_joint = gr.joint('left_thigh_joint', {-90, 0, 90}, {-90, 0, 90})
left_thigh_joint:translate(-0.75, 0.0, 0.0)
hips_joint:add_child(left_thigh_joint)

left_thigh = gr.mesh('cube', 'left_thigh')
left_thigh:scale(0.5, 1.5, 0.5)
left_thigh:translate(0.0, -1.5, 0.0)
left_thigh:set_material(red)
left_thigh_joint:add_child(left_thigh)

-- Left Calf
left_calf_joint = gr.joint('left_calf_joint', {-90, 0, 90}, {-90, 0, 90})
left_calf_joint:translate(0.0, -1.5, 0.0)
left_thigh_joint:add_child(left_calf_joint)

left_calf = gr.mesh('cube', 'left_calf')
left_calf:scale(0.5, 1.0, 0.5)
left_calf:translate(0.0, -1.0, 0.0)
left_calf:set_material(green)
left_calf_joint:add_child(left_calf)

-- Left Foot
left_foot_joint = gr.joint('left_foot_joint', {-45, 0, 45}, {-45, 0, 45})
left_foot_joint:translate(0.0, -1.0, 0.0)
left_calf_joint:add_child(left_foot_joint)

left_foot = gr.mesh('cube', 'left_foot')
left_foot:scale(0.5, 0.5, 0.5)
left_foot:translate(0.0, -0.5, 0.0)
left_foot:set_material(blue)
left_foot_joint:add_child(left_foot)

-- Right Thigh
right_thigh_joint = gr.joint('right_thigh_joint', {-90, 0, 90}, {-90, 0, 90})
right_thigh_joint:translate(0.75, 0.0, 0.0)
hips_joint:add_child(right_thigh_joint)

right_thigh = gr.mesh('cube', 'right_thigh')
right_thigh:scale(0.5, 1.5, 0.5)
right_thigh:translate(0.0, -1.5, 0.0)
right_thigh:set_material(red)
right_thigh_joint:add_child(right_thigh)

-- Right Calf
right_calf_joint = gr.joint('right_calf_joint', {-90, 0, 90}, {-90, 0, 90})
right_calf_joint:translate(0.0, -1.5, 0.0)
right_thigh_joint:add_child(right_calf_joint)

right_calf = gr.mesh('cube', 'right_calf')
right_calf:scale(0.5, 1.0, 0.5)
right_calf:translate(0.0, -1.0, 0.0)
right_calf:set_material(green)
right_calf_joint:add_child(right_calf)

-- Right Foot
right_foot_joint = gr.joint('right_foot_joint', {-45, 0, 45}, {-45, 0, 45})
right_foot_joint:translate(0.0, -1.0, 0.0)
right_calf_joint:add_child(right_foot_joint)

right_foot = gr.mesh('cube', 'right_foot')
right_foot:scale(0.5, 0.5, 0.5)
right_foot:translate(0.0, -0.5, 0.0)
right_foot:set_material(blue)
right_foot_joint:add_child(right_foot)

-- Neck
neck_joint = gr.joint('neck_joint', {-30, 0, 30}, {-30, 0, 30})
shoulder_joint:add_child(neck_joint)

neck = gr.mesh('cube', 'neck')
neck:scale(0.3, 0.5, 0.3)
neck:translate(0.0, 0.5, 0.0)
neck:set_material(white)
neck_joint:add_child(neck)

-- Head (Suzanne)
head_joint = gr.joint('head_joint', {-45, 0, 45}, {-45, 0, 45})
head_joint:translate(0.0, 0.5, 0.0)
neck_joint:add_child(head_joint)

head = gr.mesh('suzanne', 'head')
head:scale(0.8, 0.8, 0.8)
head:translate(0.0, 0.8, 0.0)  -- Adjust pivot for suzanne
head:set_material(green)
head_joint:add_child(head)

-- Move the entire root node to the correct position
root:translate(0.0, 0.0, -10.0)

return root
