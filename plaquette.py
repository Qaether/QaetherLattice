import numpy as np
import matplotlib.pyplot as plt
from itertools import permutations

def apply_d4_operations(config):
    """
    Apply all 8 D4 operations to a configuration [a,b,c,d].
    Returns a list of 8 configurations.
    
    D4 operations on a square plaquette:
    - e: identity
    - r90: rotate 90° (0->1->2->3->0)
    - r180: rotate 180° (0->2, 1->3)
    - r270: rotate 270° (0->3->2->1->0)
    - sx: reflect across horizontal (0->3, 1->2)
    - sy: reflect across vertical (0->1, 2->3)
    - sd: reflect across main diagonal (0->0, 1->3, 2->2, 3->1)
    - sd': reflect across anti-diagonal (0->2, 1->1, 2->0, 3->3)
    """
    a, b, c, d = config
    return [
        (a, b, c, d),      # e: identity
        (d, a, b, c),      # r90: rotate 90°
        (c, d, a, b),      # r180: rotate 180°
        (b, c, d, a),      # r270: rotate 270°
        (d, c, b, a),      # sx: horizontal reflection
        (b, a, d, c),      # sy: vertical reflection
        (a, d, c, b),      # sd: main diagonal reflection
        (c, b, a, d),      # sd': anti-diagonal reflection
    ]

def canonical_form(config):
    """
    Find the canonical (lexicographically smallest) form of a configuration under D4.
    """
    return min(apply_d4_operations(config))

def find_d4_equivalence_classes(values):
    """
    Find the 3 distinct equivalence classes of plaquette configurations under D4 symmetry.
    
    Args:
        values: List of 4 distinct integers
    
    Returns:
        List of 3 tuples, each representing a distinct equivalence class
    """
    # Generate all permutations
    all_perms = list(permutations(values))
    
    # Find canonical forms
    canonical_forms = set()
    for perm in all_perms:
        canonical_forms.add(canonical_form(perm))
    
    # Sort for consistent output
    return sorted(list(canonical_forms))

def generate_orthogonal_plaquettes():
    """
    Generates the coordinates for 3 orthogonal square plaquettes meeting at the origin.
    
    Returns:
        tuple: (plaquettes_corners, plaquettes_centers)
        - plaquettes_corners: List of 3 numpy arrays, each (4, 3) containing corner coordinates of a plaquette.
        - plaquettes_centers: Numpy array (3, 3) containing the center coordinates of the 3 plaquettes.
    """
    # 1. XY Plaquette (Diamond shape in XY plane)
    # Corners: (1,0,0), (0,1,0), (-1,0,0), (0,-1,0)
    corners_xy = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [-1, 0, 0],
        [0, -1, 0]
    ])
    center_xy = np.array([0, 0, 0])

    # 2. YZ Plaquette (Diamond shape in YZ plane)
    # Corners: (0,1,0), (0,0,1), (0,-1,0), (0,0,-1)
    corners_yz = np.array([
        [0, 1, 0],
        [0, 0, 1],
        [0, -1, 0],
        [0, 0, -1]
    ])
    center_yz = np.array([0, 0, 0])

    # 3. ZX Plaquette (Diamond shape in ZX plane)
    # Corners: (0,0,1), (1,0,0), (0,0,-1), (-1,0,0)
    corners_zx = np.array([
        [0, 0, 1],
        [1, 0, 0],
        [0, 0, -1],
        [-1, 0, 0]
    ])
    center_zx = np.array([0, 0, 0])

    plaquettes_corners = [corners_xy, corners_yz, corners_zx]
    plaquettes_centers = np.array([center_xy, center_yz, center_zx])

    return plaquettes_corners, plaquettes_centers

def add_edge_labels(ax, plaquettes_corners, edge_labels, colors, fontsize=14):
    """
    Add text labels on the edges of the 3 plaquettes.
    
    Args:
        ax: 3D axis
        plaquettes_corners: List of 3 corner arrays
        edge_labels: List of 3 configurations (tuples of 4 values)
        colors: List of colors for each plaquette
        fontsize: Font size for the labels (default: 14)
    """
    # Clear previous labels and arrows
    for text in ax.texts[:]:
        text.remove()
    
    # Also clear previous arrows (patches)
    for patch in ax.patches[:]:
        patch.remove()
    
    for i, (corners, config) in enumerate(zip(plaquettes_corners, edge_labels)):
        # Calculate midpoints of the 4 edges
        for j in range(4):
            start = corners[j]
            end = corners[(j + 1) % 4]
            midpoint = (start + end) / 2
            value = config[j]
            
            # Determine arrow direction based on sign
            # Positive: clockwise (from start to end)
            # Negative: counter-clockwise (from end to start)
            if value > 0:
                arrow_start = start + (end - start) * 0.25
                arrow_direction = (end - start) * 0.4
            elif value < 0:
                arrow_start = end + (start - end) * 0.25
                arrow_direction = (start - end) * 0.4
            else:  # value == 0, no arrow
                arrow_start = None
                arrow_direction = None
            
            # Add arrow if value is non-zero using quiver
            if arrow_start is not None:
                ax.quiver(arrow_start[0], arrow_start[1], arrow_start[2],
                         arrow_direction[0], arrow_direction[1], arrow_direction[2],
                         color=colors[i], arrow_length_ratio=0.3, linewidth=2)
            
            
            # Add text label
            ax.text(midpoint[0], midpoint[1], midpoint[2], 
                   str(value), 
                   fontsize=fontsize, 
                   fontweight='bold',
                   color=colors[i],
                   ha='center', 
                   va='center',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=colors[i], alpha=0.8))

def plot_orthogonal_plaquettes(plaquettes_corners, plaquettes_centers, edge_labels=None):
    """
    Plots the 3 orthogonal plaquettes.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')

    colors = ['red', 'green', 'blue']
    labels = ['XY Plane', 'YZ Plane', 'ZX Plane']

    # Plot each plaquette
    for i, corners in enumerate(plaquettes_corners):
        # Close the loop for plotting lines
        corners_plot = np.vstack([corners, corners[0]])
        
        # Plot edges
        ax.plot(corners_plot[:, 0], corners_plot[:, 1], corners_plot[:, 2], 
                color=colors[i], linewidth=2, label=f'{labels[i]} Edges')
        
        # Plot corners
        ax.scatter(corners[:, 0], corners[:, 1], corners[:, 2], 
                   color=colors[i], s=50, marker='o', alpha=0.6)

    # Plot face centers
    ax.scatter(plaquettes_centers[:, 0], plaquettes_centers[:, 1], plaquettes_centers[:, 2],
               color='black', s=100, marker='x', label='Face Centers')

    # Add origin marker
    ax.scatter([0], [0], [0], color='black', s=100, marker='o', label='Origin')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3 Orthogonal Plaquettes Forming Octahedron')
    
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_zlim([-1.2, 1.2])
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True)
    
    # Add edge labels if provided
    if edge_labels is not None:
        add_edge_labels(ax, plaquettes_corners, edge_labels, colors)

    # Adjust plot to make room for widgets
    plt.subplots_adjust(bottom=0.25)

    # Add TextBoxes
    from matplotlib.widgets import TextBox, Button, Slider

    # Status text on the plot
    status_text = plt.figtext(0.5, 0.02, "", ha="center", fontsize=12, color="red")

    # Reduced width for text boxes
    axbox1 = plt.axes([0.15, 0.15, 0.05, 0.05])
    text_box1 = TextBox(axbox1, 'Input 1: ', initial="")

    axbox2 = plt.axes([0.30, 0.15, 0.05, 0.05])
    text_box2 = TextBox(axbox2, 'Input 2: ', initial="")

    axbox3 = plt.axes([0.45, 0.15, 0.05, 0.05])
    text_box3 = TextBox(axbox3, 'Input 3: ', initial="")

    axbox4 = plt.axes([0.60, 0.15, 0.05, 0.05])
    text_box4 = TextBox(axbox4, 'Input 4: ', initial="")

    # Add Font Size Slider
    axslider = plt.axes([0.15, 0.08, 0.4, 0.03])
    font_slider = Slider(axslider, 'Font Size', 6, 24, valinit=14, valstep=1)
    
    # Store current equivalence classes
    current_classes = [None]
    
    def update_fontsize(val):
        if current_classes[0] is not None:
            add_edge_labels(ax, plaquettes_corners, current_classes[0], colors, fontsize=int(val))
            plt.draw()
    
    font_slider.on_changed(update_fontsize)

    # Add Run Button
    axrun = plt.axes([0.8, 0.05, 0.1, 0.075])
    b_run = Button(axrun, 'Run')

    def run_callback(event):
        inputs = [text_box1.text, text_box2.text, text_box3.text, text_box4.text]
        try:
            # 1. Parse integers
            values = [int(val) for val in inputs]
            
            # 2. Check range [-5, 6]
            for v in values:
                if not (-5 <= v <= 6):
                    status_text.set_text(f"Error: Value {v} out of range [-5, 6].")
                    status_text.set_color("red")
                    plt.draw()
                    return

            # 3. Check sum {-12, 0, 12}
            total_sum = sum(values)
            if total_sum not in [-12, 0, 12]:
                status_text.set_text(f"Error: Sum is {total_sum}. Must be -12, 0, or 12.")
                status_text.set_color("red")
                plt.draw()
                return

            # 4. Find 3 distinct D4 equivalence classes
            equivalence_classes = find_d4_equivalence_classes(values)
            
            # Store for slider updates
            current_classes[0] = equivalence_classes
            
            # 5. Add edge labels to the plaquettes with current fontsize
            current_fontsize = int(font_slider.val)
            add_edge_labels(ax, plaquettes_corners, equivalence_classes, colors, fontsize=current_fontsize)
            
            # Format output
            config_strs = [f"[{','.join(map(str, config))}]" for config in equivalence_classes]
            result_text = f"Success! 3 distinct D4 configurations: {' | '.join(config_strs)}"
            
            status_text.set_text(result_text)
            status_text.set_color("green")
            print(f"Success! Valid inputs: {values}, Sum: {total_sum}")
            print(f"3 D4 equivalence classes: {equivalence_classes}")
            plt.draw()
            
        except ValueError:
            status_text.set_text("Error: All inputs must be valid integers.")
            status_text.set_color("red")
            plt.draw()

    b_run.on_clicked(run_callback)
    
    plt.show()

if __name__ == "__main__":
    corners, centers = generate_orthogonal_plaquettes()
    plot_orthogonal_plaquettes(corners, centers)
