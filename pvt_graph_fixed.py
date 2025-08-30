import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import streamlit as st
import matplotlib.animation as animation
from io import BytesIO
import base64

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


def create_temperature_animation(num_frames, a_val, b_val, n_val, p_max_val):
    """Create an animated GIF of temperature sweep"""
    # Temperature range for animation
    temps = np.linspace(200, 400, num_frames)
    
    # Volume range
    V_range = np.linspace(0.00001, 1, 100)  # Reduced resolution for performance
    
    # Set up the figure and axis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Initialize empty plots
    line1, = ax1.plot([], [], 'red', linewidth=3)
    line2, = ax2.plot([], [], 'red', linewidth=3)
    
    def animate_frame(frame):
        T_current = temps[frame]
        
        # Calculate pressure for current temperature
        P_slice = van_der_waals_equation(V_range, T_current, a_val, b_val, n_val)
        P_slice = np.where(P_slice >= 0, P_slice, np.nan)
        
        # Update left plot (P-V)
        valid_mask = ~np.isnan(P_slice)
        if np.any(valid_mask):
            line1.set_data(V_range[valid_mask], P_slice[valid_mask])
        else:
            line1.set_data([], [])
        
        # Update right plot (same as left for now, could be different view)
        if np.any(valid_mask):
            line2.set_data(V_range[valid_mask], P_slice[valid_mask])
        else:
            line2.set_data([], [])
        
        # Update titles
        ax1.set_title(f'P-V Isotherme bei T = {T_current:.0f} K')
        ax2.set_title(f'P-V Isotherme (Zoom) bei T = {T_current:.0f} K')
        
        return line1, line2
    
    # Set up axes
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, p_max_val)
    ax1.set_xlabel('Volumen (V) [L]')
    ax1.set_ylabel('Druck (P) [Pa]')
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlim(0, 0.2)  # Zoomed view
    ax2.set_ylim(0, p_max_val * 0.5)
    ax2.set_xlabel('Volumen (V) [L]')
    ax2.set_ylabel('Druck (P) [Pa]')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate_frame, frames=num_frames, 
                                 interval=200, blit=True, repeat=True)
    
    # Convert to HTML
    html_str = anim.to_jshtml()
    plt.close(fig)
    
    return html_str
    

def create_3d_pbt_diagram():
    # Streamlit page configuration
    st.set_page_config(page_title="Van der Waals 3D PVT Diagramm", layout="wide")
    st.title("3D PVT Diagramm - Van der Waals Zustandsgleichung")
    
    # Sidebar for controls
    st.sidebar.header("Parameter")
    
    # Initial parameters for Van der Waals equation - CO₂ values
    T_init = 300
    a_init = 364.0   # Van der Waals parameter a for CO₂ [Pa·L²/mol²]
    b_init = 0.04267  # Van der Waals parameter b for CO₂ [L/mol]
    n_init = 1.0  # Changed to float
    P_max_init = 50000
    
    # Streamlit controls in sidebar
    T_slice_val = st.sidebar.slider('Temperatur-Schnitt (K)', min_value=200, max_value=400, value=T_init, step=1)
    
    # Text inputs for a and b parameters
    col1, col2 = st.sidebar.columns(2)
    with col1:
        a_val = st.number_input('Parameter a', value=a_init, format="%.1f", key="param_a")
    with col2:
        b_val = st.number_input('Parameter b', value=b_init, format="%.5f", key="param_b")
    
    n_val = st.sidebar.slider('Stoffmenge n (mol)', min_value=0.1, max_value=5.0, value=n_init, step=0.1)
    p_max_val = st.sidebar.slider('P-Achse Max (Pa)', min_value=1000, max_value=100000, value=P_max_init, step=1000)
    
    # Animation controls
    st.sidebar.markdown("---")
    st.sidebar.subheader("Animation")
    
    # Animation button and speed control
    col_anim1, col_anim2 = st.sidebar.columns(2)
    with col_anim1:
        generate_animation = st.button("🎬 Animation generieren", key="generate_anim")
    with col_anim2:
        animation_frames = st.selectbox("Frames", [20, 40, 60, 80], index=1, key="frames")
    
    # Static animation display
    if generate_animation or 'animation_html' in st.session_state:
        if generate_animation:
            # Generate animation only when button is clicked
            with st.spinner('Generiere Animation... Bitte warten.'):
                animation_html = create_temperature_animation(animation_frames, a_val, b_val, n_val, p_max_val)
                st.session_state.animation_html = animation_html
        
        if 'animation_html' in st.session_state:
            st.sidebar.success("✅ Animation bereit!")
            st.sidebar.markdown("**Animation wird unten angezeigt**")
    
    # Create volume and temperature ranges with high resolution
    V_range = np.linspace(0.00001, 1, 200)
    T_range = np.linspace(200, 400, 200)
    V_mesh, T_mesh = np.meshgrid(V_range, T_range)

    # Calculate pressure surface
    P_mesh = van_der_waals_equation(V_mesh, T_mesh, a_val, b_val, n_val)
    P_mesh = np.where(P_mesh >= 0, P_mesh, np.nan)  # Mask negative pressures

    # Calculate pressure slice for current temperature
    V_slice_range = V_range
    P_slice = van_der_waals_equation(V_slice_range, T_slice_val, a_val, b_val, n_val)
    P_slice = np.where(P_slice >= 0, P_slice, np.nan)  # Mask negative pressures

    # Create two columns for the plots
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("3D PVT Oberfläche")
        
        # Create 3D plot
        fig_3d = plt.figure(figsize=(10, 8))
        ax3d = fig_3d.add_subplot(111, projection='3d')
        
        # Create 3D surface plot
        surface = ax3d.plot_surface(
            V_mesh, T_mesh, P_mesh, cmap='coolwarm', alpha=0.7,
            linewidth=0, antialiased=False
        )
        
        # Add temperature slice line - ensure it's always visible on top
        valid_mask = ~np.isnan(P_slice)
        if np.any(valid_mask):
            ax3d.plot(V_slice_range[valid_mask],
                     np.full_like(V_slice_range[valid_mask], T_slice_val),
                     P_slice[valid_mask], 'red', linewidth=6, 
                     zorder=1000, alpha=1.0)
        
        # Set labels and limits
        ax3d.set_xlabel('Volumen (V) [L]')
        ax3d.set_ylabel('Temperatur (T) [K]')
        ax3d.set_zlabel('Druck (P) [Pa]')
        ax3d.set_zlim(0, p_max_val)
        
        # Add animation indicator to title
        ax3d.set_title(f'Van der Waals: a={a_val:.0f}, b={b_val:.3f}, n={n_val:.1f}')
        
        st.pyplot(fig_3d)
    
    with col2:
        st.subheader("P-V Isotherme (2D Schnitt)")
        
        # Create 2D plot
        fig_2d = plt.figure(figsize=(10, 8))
        ax2d = fig_2d.add_subplot(111)
        
        # Plot 2D slice
        valid_mask_2d = ~np.isnan(P_slice)
        if np.any(valid_mask_2d):
            ax2d.plot(V_slice_range[valid_mask_2d], P_slice[valid_mask_2d], 'red', linewidth=3)
        
        ax2d.set_xlabel('Volumen (V) [L]')
        ax2d.set_ylabel('Druck (P) [Pa]')
        ax2d.set_ylim(0, p_max_val)
        ax2d.set_title(f'P-V Isotherme bei T = {T_slice_val:.0f} K')
        ax2d.grid(True, alpha=0.3)
        
        st.pyplot(fig_2d)
    
    # Display animation if available
    if 'animation_html' in st.session_state:
        st.markdown("---")
        st.subheader("🎬 Temperatur-Animation")
        st.components.v1.html(st.session_state.animation_html, height=500, scrolling=True)
    
    # Add information section
    st.markdown("---")
    st.subheader("Anleitung")
    st.markdown("""
    - **Temperatur-Schnitt**: Verwende den Schieberegler, um die Temperatur der roten Linie im 3D-Plot zu ändern
    - **Parameter a & b**: Gib präzise Van der Waals Parameter für verschiedene Gase ein
    - **Stoffmenge n**: Passe die Anzahl der Mol an
    - **P-Achse Max**: Ändere den maximalen Druck für bessere Visualisierung
    - **🎬 Animation generieren**: Erstellt eine flüssige HTML5-Animation (einmalig generieren, dann wiederverwenden)
    
    **Animations-Features:**
    - **Effizient für Server**: Keine ständigen Neuladezyklen
    - **Einstellbare Frames**: 20-80 Frames für verschiedene Qualitätsstufen
    - **Flüssige Wiedergabe**: HTML5-basierte Animation läuft im Browser
    - **Wiederverwendbar**: Einmal generiert, kann beliebig oft abgespielt werden
    - **Zwei Ansichten**: Vollansicht und Zoom-Ansicht der P-V Isotherme
    
    **Aktuelle CO₂ Parameter:**
    - a = 364,0 Pa·L²/mol² (zwischenmolekulare Anziehung)
    - b = 0,04267 L/mol (Molekülvolumen)
    """)
    
    # Close matplotlib figures to prevent memory leaks
    plt.close('all')

if __name__ == "__main__":
    create_3d_pbt_diagram()
