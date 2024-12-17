import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D
import math
import time

class GeometryCalculator:
    def __init__(self, root):
        self.root = root
        self.root.title("Geometry Calculator")
        
        # main container
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # left side - controls
        self.control_frame = ttk.Frame(self.main_frame)
        self.control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5)
        
        # list of shapes
        self.shapes = ["Cube", "Cuboid", "Cylinder", "Cone", "Sphere"]
        self.shape_var = tk.StringVar()
        self.shape_combo = ttk.Combobox(self.control_frame, 
                                      textvariable=self.shape_var,
                                      values=self.shapes)
        self.shape_combo.pack(pady=5)
        self.shape_combo.bind('<<ComboboxSelected>>', self.on_shape_selected)
        
        # container for the dimensions fields
        self.dimensions_frame = ttk.Frame(self.control_frame)
        self.dimensions_frame.pack(fill=tk.X, pady=5)
        
        # right side - visualization
        self.fig = plt.Figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_frame)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        # dictionary storing the fields for the dimensions
        self.dimension_entries = {}
        
        # dictionary storing the definitions of the dimensions for each shape
        self.shape_dimensions = {
            "Cube": [("Side", "a")],
            "Cuboid": [
                ("Length", "a"),
                ("Width", "b"),
                ("Height", "c")
            ],
            "Cylinder": [
                ("Radius", "r"),
                ("Height", "h")
            ],
            "Cone": [
                ("Radius", "r"),
                ("Height", "h")
            ],
            "Sphere": [("Radius", "r")]
        }
        
        # variables for rotation
        self.rotation_active = True
        self.last_interaction_time = 0
        self.rotation_pause_duration = 3  # seconds
        
        # adding mouse handling
        self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
        self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        
        self.is_rotating = False
        self.prev_x = 0
        self.prev_y = 0

    def on_shape_selected(self, event=None):
        # clearing the previous fields for the dimensions
        for widget in self.dimensions_frame.winfo_children():
            widget.destroy()
        self.dimension_entries.clear()
        
        # creating new fields for the selected shape
        shape = self.shape_var.get()
        if shape in self.shape_dimensions:
            for label, key in self.shape_dimensions[shape]:
                self.create_dimension_entry(label, key)
        
        # force the redraw of the visualization
        self.ax.clear()
        self.set_plot_properties()

    def create_dimension_entry(self, label, key):
        frame = ttk.Frame(self.dimensions_frame)
        frame.pack(fill=tk.X, pady=2)
        
        ttk.Label(frame, text=label).pack(side=tk.LEFT)
        entry = ttk.Entry(frame, width=10)
        entry.pack(side=tk.LEFT, padx=5)
        entry.bind('<KeyRelease>', self.update_visualization)
        
        self.dimension_entries[key] = entry 

    def draw_cube(self, side, scale):
        # scale the dimensions while keeping the proportions
        a = side * scale
        
        # defining the vertices
        points = np.array([
            [0, 0, 0], [a, 0, 0], [a, a, 0], [0, a, 0],
            [0, 0, a], [a, 0, a], [a, a, a], [0, a, a]
        ])
        
        edges = [[0,1], [1,2], [2,3], [3,0],
                [4,5], [5,6], [6,7], [7,4],
                [0,4], [1,5], [2,6], [3,7]]
                
        for edge in edges:
            self.ax.plot3D(*zip(*points[edge]), color='blue')
            
        # add the label for the dimension
        self.ax.text(a/2, -0.1, 0, f'a={side:.1f}', ha='center')

    def start_rotation(self):
        self.rotation_active = True
        self.rotate()

    def rotate(self):
        current_time = time.time()
        
        # check if 10 seconds have passed since the last interaction
        if not self.rotation_active and \
           (current_time - self.last_interaction_time) > self.rotation_pause_duration:
            self.rotation_active = True
        
        if self.rotation_active:
            self.ax.azim = (self.ax.azim + 2) % 360
            self.canvas.draw()
            
        self.root.after(50, self.rotate)  # continue checking every 50ms

    def update_visualization(self, event=None):
        try:
            shape = self.shape_var.get()
            dimensions = self.get_dimensions()
            if not dimensions:
                return
                
            # constant maximum length for the largest dimension
            MAX_DISPLAY_SIZE = 1.0
            
            # find the largest dimension to scale
            max_dim = max(dimensions.values())
            scale_factor = MAX_DISPLAY_SIZE / max_dim if max_dim > 0 else 1.0
            
            self.ax.clear()
            
            # use only the dimensions that are needed for the current shape
            if shape == "Cube" and 'a' in dimensions:
                self.draw_cube(dimensions['a'], scale_factor)
            elif shape == "Cuboid" and all(k in dimensions for k in ['a', 'b', 'c']):
                self.draw_cuboid(dimensions['a'], dimensions['b'], dimensions['c'], scale_factor)
            elif shape == "Cylinder" and all(k in dimensions for k in ['r', 'h']):
                self.draw_cylinder(dimensions['r'], dimensions['h'], scale_factor)
            elif shape == "Cone" and all(k in dimensions for k in ['r', 'h']):
                self.draw_cone(dimensions['r'], dimensions['h'], scale_factor)
            elif shape == "Sphere" and 'r' in dimensions:
                self.draw_sphere(dimensions['r'], scale_factor)
                
            self.set_plot_properties()
            self.calculate_and_display_metrics(dimensions)
            
        except (ValueError, KeyError):
            pass

    def get_dimensions(self):
        dimensions = {}
        shape = self.shape_var.get()
        # get only the dimensions that are defined for the current shape
        required_dims = [dim[1] for dim in self.shape_dimensions[shape]]
        
        for key, entry in self.dimension_entries.items():
            if key in required_dims:  # check if the dimension is needed for the current shape
                try:
                    value = float(entry.get())
                    if value <= 0:
                        return None
                    dimensions[key] = value
                except ValueError:
                    return None
        return dimensions if len(dimensions) == len(required_dims) else None

    def draw_cuboid(self, length, width, height, scale):
        # scale the dimensions while keeping the proportions
        l = length * scale
        w = width * scale
        h = height * scale
        
        # define the vertices
        points = np.array([
            [0, 0, 0], [l, 0, 0], [l, w, 0], [0, w, 0],
            [0, 0, h], [l, 0, h], [l, w, h], [0, w, h]
        ])
        
        edges = [[0,1], [1,2], [2,3], [3,0],
                 [4,5], [5,6], [6,7], [7,4],
                 [0,4], [1,5], [2,6], [3,7]]
                 
        for edge in edges:
            self.ax.plot3D(*zip(*points[edge]), color='blue')
        
        # add the labels for the dimensions
        self.ax.text(l/2, -0.1, 0, f'{length:.1f}', ha='center')
        self.ax.text(l+0.1, w/2, 0, f'{width:.1f}', ha='left')
        self.ax.text(0, 0, h/2, f'{height:.1f}', va='center')

    def draw_cylinder(self, radius, height, scale):
        r = radius * scale
        h = height * scale
        
        # creating points for the base circles
        theta = np.linspace(0, 2*np.pi, 100)
        
        # bottom base
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        z = np.zeros_like(theta)
        self.ax.plot(x, y, z, 'blue')
        
        # top base
        self.ax.plot(x, y, np.full_like(theta, h), 'blue')
        
        # lines connecting the bases
        for angle in np.linspace(0, 2*np.pi, 16):
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            self.ax.plot([x,x], [y,y], [0,h], 'blue')
        
        # labels for the dimensions
        self.ax.text(r+0.1, 0, 0, f'r={radius:.1f}')
        self.ax.text(0, 0, h/2, f'h={height:.1f}')

    def draw_cone(self, radius, height, scale):
        r = radius * scale
        h = height * scale
        
        # base
        theta = np.linspace(0, 2*np.pi, 100)
        x_base = r * np.cos(theta)
        y_base = r * np.sin(theta)
        z_base = np.zeros_like(theta)
        
        self.ax.plot(x_base, y_base, z_base, 'blue')
        
        # lines to the vertex
        for angle in np.linspace(0, 2*np.pi, 16):
            x = r * np.cos(angle)
            y = r * np.sin(angle)
            self.ax.plot([x,0], [y,0], [0,h], 'blue')
        
        # labels for the dimensions
        self.ax.text(r+0.1, 0, 0, f'r={radius:.1f}')
        self.ax.text(0, 0, h/2, f'h={height:.1f}')

    def draw_sphere(self, radius, scale):
        r = radius * scale
        
        # creating a grid of points for the sphere
        phi = np.linspace(0, np.pi, 20)
        theta = np.linspace(0, 2*np.pi, 40)
        phi, theta = np.meshgrid(phi, theta)
        
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        
        self.ax.plot_wireframe(x, y, z, color='blue', alpha=0.5)
        
        # labels for the dimensions
        self.ax.text(r+0.1, 0, 0, f'r={radius:.1f}')

    def set_plot_properties(self):
        # set the equal proportions for the axes
        self.ax.set_box_aspect([1,1,1])
        
        # find the maximum range for all axes
        all_limits = [self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()]
        max_range = max(abs(max(x)) for x in all_limits)
        
        # set the symmetrical limits for all axes
        self.ax.set_xlim(-max_range, max_range)
        self.ax.set_ylim(-max_range, max_range)
        self.ax.set_zlim(-max_range, max_range)
        
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')
        self.ax.grid(True)
        self.canvas.draw()

    def calculate_and_display_metrics(self, dimensions):
        shape = self.shape_var.get()
        volume = 0
        surface_area = 0
        
        if shape == "Cube":
            a = dimensions['a']
            volume = a**3
            surface_area = 6 * a**2
        elif shape == "Cuboid":
            a, b, c = dimensions['a'], dimensions['b'], dimensions['c']
            volume = a * b * c
            surface_area = 2 * (a*b + b*c + a*c)
        elif shape == "Cylinder":
            r, h = dimensions['r'], dimensions['h']
            volume = np.pi * r**2 * h
            surface_area = 2 * np.pi * r**2 + 2 * np.pi * r * h
        elif shape == "Cone":
            r, h = dimensions['r'], dimensions['h']
            volume = (1/3) * np.pi * r**2 * h
            l = np.sqrt(r**2 + h**2)  # tworząca
            surface_area = np.pi * r**2 + np.pi * r * l
        elif shape == "Sphere":
            r = dimensions['r']
            volume = (4/3) * np.pi * r**3
            surface_area = 4 * np.pi * r**2
            
        # update the labels with the results
        if not hasattr(self, 'results_frame'):
            self.results_frame = ttk.Frame(self.control_frame)
            self.results_frame.pack(fill=tk.X, pady=10)
            self.volume_label = ttk.Label(self.results_frame, text="")
            self.volume_label.pack()
            self.area_label = ttk.Label(self.results_frame, text="")
            self.area_label.pack()
            
        self.volume_label.config(text=f"Volume: {volume:.2f}")
        self.area_label.config(text=f"Surface Area: {surface_area:.2f}")

    def on_mouse_press(self, event):
        if event.inaxes == self.ax:
            self.is_rotating = True
            self.prev_x = event.xdata
            self.prev_y = event.ydata
            self.pause_auto_rotation()

    def on_mouse_release(self, event):
        self.is_rotating = False
        self.last_interaction_time = time.time()

    def on_mouse_move(self, event):
        if self.is_rotating and event.inaxes == self.ax:
            dx = event.xdata - self.prev_x
            dy = event.ydata - self.prev_y
            
            # updating the view angles
            current_elev = self.ax.elev
            current_azim = self.ax.azim
            
            # changing the view angle based on the mouse movement
            self.ax.view_init(elev=current_elev + dy * 180,
                            azim=current_azim - dx * 180)
            
            self.prev_x = event.xdata
            self.prev_y = event.ydata
            self.canvas.draw()

    def pause_auto_rotation(self):
        self.rotation_active = False
        self.last_interaction_time = time.time()

if __name__ == "__main__":
    root = tk.Tk()
    app = GeometryCalculator(root)
    app.start_rotation()  #start the rotation   
    root.mainloop()