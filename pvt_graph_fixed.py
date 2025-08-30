import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, TextBox
import time

def van_der_waals_equation(V, T, a, b, n, p_max=None):
    """
    Van der Waals equation of state
    P = (nRT)/(V - nb) - (an²)/V²
    
    Units:
    - V: volume [L]
    - T: temperature [K] 
    - a: [Pa·L²/mol²]
    - b: [L/mol]
    - n: moles [mol]
    - p_max: maximum pressure to display [Pa] (optional)
    - P: pressure [Pa]
    """
    R = 8.314  # Universal gas constant [J/(mol·K)]
    
    # Avoid division by zero
    V_safe = np.maximum(V, 0.00001)
    T_safe = np.maximum(T, 1.0)
    
    # Check physical constraints (V > nb)
    condition = V_safe > n * b
    
    # Calculate Van der Waals pressure
    P_rep = (n * R * T_safe) / (V_safe - n * b)
    P_att = (a * n**2) / (V_safe**2)
    
    P_calc = P_rep - P_att
    
    # Only return pressures where conditions are met
    P = np.where(condition, P_calc, np.nan)  # Use NaN for invalid regions
    
 
    
    return P
    

def create_3d_pbt_diagram():
    # Initial parameters for Van der Waals equation - CO₂ values
    T_init = 300
    a_init = 364.0   # Van der Waals parameter a for CO₂ [Pa·L²/mol²]
    b_init = 0.04267  # Van der Waals parameter b for CO₂ [L/mol]
    n_init = 1

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    ax3d = fig.add_subplot(121, projection='3d')
    ax2d = fig.add_subplot(122)
    plt.subplots_adjust(left=0.05, bottom=0.25, right=0.95, top=0.9, wspace=0.3)

    # Create volume and temperature ranges with high resolution
    V_range = np.linspace(0.00001, 1, 200)
    T_range = np.linspace(200, 400, 200)
    V_mesh, T_mesh = np.meshgrid(V_range, T_range)

    # Calculate initial pressure surface
    P_mesh = van_der_waals_equation(V_mesh, T_mesh, a_init, b_init, n_init)
    P_mesh = np.where(P_mesh >= 0, P_mesh, np.nan)  # Mask negative pressures

    # Create 3D surface plot without grid lines
    surface = ax3d.plot_surface(
        V_mesh, T_mesh, P_mesh, cmap='coolwarm', alpha=0.7,
        linewidth=0, antialiased=False
    )

    # Add temperature slice - ensure it's always visible on top
    T_slice = T_init
    V_slice_range = V_range
    P_slice = van_der_waals_equation(V_slice_range, T_slice, a_init, b_init, n_init)
    P_slice = np.where(P_slice >= 0, P_slice, np.nan)  # Mask negative pressures
    slice_line = ax3d.plot(V_slice_range[~np.isnan(P_slice)],
                          np.full_like(V_slice_range[~np.isnan(P_slice)], T_slice),
                          P_slice[~np.isnan(P_slice)], 'red', linewidth=3, 
                          zorder=1000, alpha=1.0)[0]  # High zorder and full opacity

    # Set labels (no axis inversion)
    ax3d.set_xlabel('Volume (V) [L]')
    ax3d.set_ylabel('Temperature (T) [K]')
    ax3d.set_zlabel('Pressure (P) [Pa]')
    ax3d.set_title('3D PVT Diagram - Van der Waals EOS')

    # Plot 2D slice
    line2d, = ax2d.plot(V_slice_range[~np.isnan(P_slice)], P_slice[~np.isnan(P_slice)], 'red', linewidth=3)
    ax2d.set_xlabel('Volume (V) [L]')
    ax2d.set_ylabel('Pressure (P) [Pa]')
    ax2d.set_title('P-V Isotherm (Temperature Slice)')
    ax2d.grid(True, alpha=0.3)

    # Create sliders and text boxes
    ax_T_slice = plt.axes([0.1, 0.18, 0.3, 0.02])
    ax_a = plt.axes([0.1, 0.15, 0.15, 0.02])  # Smaller width for text box
    ax_b = plt.axes([0.25, 0.15, 0.15, 0.02])  # Side by side with a
    ax_n = plt.axes([0.1, 0.09, 0.3, 0.02])
    ax_p_max = plt.axes([0.1, 0.06, 0.3, 0.02])

    # Set initial pressure axis maximum to 50000
    P_max_init = 50000

    slider_T_slice = Slider(ax_T_slice, 'Slice T', 200, 400, valinit=T_init)
    textbox_a = TextBox(ax_a, 'Parameter a: ', initial=str(a_init))
    textbox_b = TextBox(ax_b, 'Parameter b: ', initial=str(b_init))
    slider_n = Slider(ax_n, 'Moles n', 0.1, 5, valinit=n_init)
    slider_p_max = Slider(ax_p_max, 'P-axis Max', 1000, 100000, valinit=P_max_init)

    # Optimization: Add simple debouncing to prevent too many rapid updates
    last_update_time = [0]  # Use list to make it mutable in nested function
    
    # Add keyboard controls for temperature slicing
    def on_key_press(event):
        """Handle keyboard input for temperature control"""
        if event.key == 'up':
            # Increase temperature by 5K
            new_temp = min(slider_T_slice.val + 5, slider_T_slice.valmax)
            slider_T_slice.set_val(new_temp)
        elif event.key == 'down':
            # Decrease temperature by 5K
            new_temp = max(slider_T_slice.val - 5, slider_T_slice.valmin)
            slider_T_slice.set_val(new_temp)
        elif event.key == 'right':
            # Increase temperature by 1K for fine control
            new_temp = min(slider_T_slice.val + 1, slider_T_slice.valmax)
            slider_T_slice.set_val(new_temp)
        elif event.key == 'left':
            # Decrease temperature by 1K for fine control
            new_temp = max(slider_T_slice.val - 1, slider_T_slice.valmin)
            slider_T_slice.set_val(new_temp)
    
    # Connect keyboard event handler
    fig.canvas.mpl_connect('key_press_event', on_key_press)

    
    def update_plot(val):
        nonlocal surface, slice_line  # Declare nonlocal at the beginning
        
        # Enhanced debouncing: prevent updates more frequent than 30ms for smoother feel
        current_time = time.time()
        if current_time - last_update_time[0] < 0.01:
            return
        last_update_time[0] = current_time
        
        T_slice_val = slider_T_slice.val
        # Get values from text boxes with error handling
        try:
            a_val = float(textbox_a.text)
        except ValueError:
            a_val = a_init  # Use default if invalid input
        try:
            b_val = float(textbox_b.text)
        except ValueError:
            b_val = b_init  # Use default if invalid input
        n_val = slider_n.val
        p_max_val = slider_p_max.val
        
        # Calculate new pressure surface more efficiently
        P_mesh_new = van_der_waals_equation(V_mesh, T_mesh, a_val, b_val, n_val, p_max_val)
        P_slice_new = van_der_waals_equation(V_slice_range, T_slice_val, a_val, b_val, n_val, p_max_val)
        
        # Update 3D surface data without clearing (much faster)
        try:
            surface.remove()
        except:
            pass  # Handle case where surface might already be removed
        surface = ax3d.plot_surface(
            V_mesh, T_mesh, P_mesh_new, cmap='coolwarm', alpha=0.7,
            linewidth=0, antialiased=False, zorder=1  # Lower z-order than the red line
        )
        
        # Update the temperature slice line efficiently - ensure it stays on top
        valid_mask = ~np.isnan(P_slice_new)
        if np.any(valid_mask):
            # Remove and recreate the line to ensure it's on top
            slice_line.remove()
            slice_line = ax3d.plot(V_slice_range[valid_mask], 
                                  np.full_like(V_slice_range[valid_mask], T_slice_val), 
                                  P_slice_new[valid_mask], 'red', linewidth=3, 
                                  zorder=1000, alpha=1.0)[0]  # Always on top
        else:
            slice_line.set_visible(False)
        
        # Update axis limits and title without full redraw
        ax3d.set_zlim(0, p_max_val)
        ax3d.set_title(f'Van der Waals: a={a_val:.0f}, b={b_val:.3f}, n={n_val:.1f}')
        
        # Update 2D plot more efficiently
        valid_mask_2d = ~np.isnan(P_slice_new)
        if np.any(valid_mask_2d):
            line2d.set_data(V_slice_range[valid_mask_2d], P_slice_new[valid_mask_2d])
        else:
            line2d.set_data([], [])
        
        # Update 2D plot limits and title
        ax2d.set_ylim(0, p_max_val)
        ax2d.set_title(f'P-V Isotherm at T = {T_slice_val:.0f} K')
        
        # Use even more efficient drawing
        fig.canvas.draw_idle()
    
    slider_T_slice.on_changed(update_plot)
    textbox_a.on_submit(update_plot)  # TextBox uses on_submit instead of on_changed
    textbox_b.on_submit(update_plot)  # TextBox uses on_submit instead of on_changed
    slider_n.on_changed(update_plot)
    slider_p_max.on_changed(update_plot)  # Add event handler for pressure axis control
    
    # Add instructions text
    fig.text(0.02, 0.02, 'Keyboard Controls:\n↑/↓ arrows: ±5K temperature\n←/→ arrows: ±1K temperature\nClick plot to activate keyboard control', 
             fontsize=9, bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    plt.show()

if __name__ == "__main__":
    print("🚀 Starting 3D PVT Diagram with Van der Waals Equation")
    print("📋 Keyboard Controls:")
    print("   ↑/↓ Arrow Keys: Change temperature by ±5K")
    print("   ←/→ Arrow Keys: Change temperature by ±1K")
    print("   💡 Tip: Click on the plot area first to activate keyboard controls!")
    create_3d_pbt_diagram()
