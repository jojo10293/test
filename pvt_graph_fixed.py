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
    """Create animated plots showing temperature sweep through Van der Waals surface"""
    # Temperature range for animation
    temps = np.linspace(200, 400, num_frames)
    
    # Create volume and temperature ranges for 3D surface
    V_range = np.linspace(0.00001, 1, 80)  # Reduced resolution for performance
    T_range = np.linspace(200, 400, 80)
    V_mesh, T_mesh = np.meshgrid(V_range, T_range)
    
    # Calculate the full 3D surface once
    P_mesh = van_der_waals_equation(V_mesh, T_mesh, a_val, b_val, n_val)
    #P_mesh = np.where(P_mesh >= 0, P_mesh, np.nan)
    
    # Set up figure with fixed, larger size
    fig = plt.figure(figsize=(16*0.8, 8*0.8))
    fig.set_tight_layout(True)
    
    # 3D subplot
    ax3d = fig.add_subplot(121, projection='3d')
    ax3d.plot_surface(V_mesh, T_mesh, P_mesh, cmap='coolwarm', alpha=0.6, 
                     linewidth=0, antialiased=False)
    line3d, = ax3d.plot([], [], [], 'red', linewidth=3, alpha=1.0, zorder=1000)
    
    ax3d.set_xlabel('Volumen (V) [L]', fontsize=10)
    ax3d.set_ylabel('Temperatur (T) [K]', fontsize=10)
    ax3d.set_zlabel('Druck (P) [Pa]', fontsize=10)
    ax3d.set_xlim(0, 1)
    ax3d.set_ylim(200, 400)
    ax3d.set_zlim(0, p_max_val)
    
    # 2D subplot
    ax2d = fig.add_subplot(122)
    line2d, = ax2d.plot([], [], 'red', linewidth=3)
    ax2d.set_xlabel('Volumen (V) [L]', fontsize=12)
    ax2d.set_ylabel('Druck (P) [Pa]', fontsize=12)
    ax2d.set_xlim(0, 1)
    ax2d.set_ylim(0, p_max_val)
    ax2d.grid(True, alpha=0.3)
    
    def animate_frame(frame):
        T_current = temps[frame]
        
        # Calculate pressure slice for current temperature
        P_slice = van_der_waals_equation(V_range, T_current, a_val, b_val, n_val)
        #P_slice = np.where(P_slice >= 0, P_slice, np.nan)
        
        valid_mask = ~np.isnan(P_slice)
        
        if np.any(valid_mask):
            # Update 3D line
            line3d.set_data_3d(V_range[valid_mask], 
                              np.full_like(V_range[valid_mask], T_current), 
                              P_slice[valid_mask])
            
            # Update 2D line
            line2d.set_data(V_range[valid_mask], P_slice[valid_mask])
        else:
            line3d.set_data_3d([], [], [])
            line2d.set_data([], [])
        
        # Update titles
        ax3d.set_title(f'Van der Waals 3D\nT = {T_current:.0f} K', fontsize=12)
        ax2d.set_title(f'P-V Isotherme\nT = {T_current:.0f} K', fontsize=12)
        
        return line3d, line2d
    
    # Create animation with better settings for stability
    anim = animation.FuncAnimation(fig, animate_frame, frames=num_frames, 
                                 interval=400, blit=False, repeat=True, cache_frame_data=False)
    
    # Convert to HTML with responsive CSS styling
    html_str = anim.to_jshtml(fps=2.5, embed_frames=True, default_mode='loop')
    
    # Simple CSS without responsive complications
    responsive_html = f"""
    <style>
        .animation-container {{
            width: 100%;
            text-align: center;
        }}
    </style>
    <div class="animation-container">
        {html_str}
    </div>
    """
    
    plt.close(fig)
    
    return responsive_html
    

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
    T_slice_val = st.sidebar.slider('Temperatur-Schnitt (K)', min_value=200, max_value=400, value=T_init, step=1, key="temp_slider")
    
    # Text inputs for a and b parameters
    col1, col2 = st.sidebar.columns(2)
    with col1:
        a_val = st.number_input('Parameter a', value=a_init, format="%.1f", key="param_a_input")
    with col2:
        b_val = st.number_input('Parameter b', value=b_init, format="%.5f", key="param_b_input")
    
    n_val = st.sidebar.slider('Stoffmenge n (mol)', min_value=0.1, max_value=5.0, value=n_init, step=0.1, key="n_slider")
    p_max_val = st.sidebar.slider('P-Achse Max (Pa)', min_value=1000, max_value=100000, value=P_max_init, step=1000, key="p_max_slider")
    
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
    
    # Create cache key for surface data
    cache_key = f"{a_val}_{b_val}_{n_val}_{p_max_val}"
    
    # Cache expensive calculations
    if 'surface_cache_key' not in st.session_state or st.session_state.surface_cache_key != cache_key:
        # Only recalculate if parameters changed
        with st.spinner('Berechne Oberfläche...'):
            # Create volume and temperature ranges with optimized resolution
            V_range = np.linspace(0.00001, 1, 150)  # Reduced from 200 for speed
            T_range = np.linspace(200, 400, 150)
            V_mesh, T_mesh = np.meshgrid(V_range, T_range)

            # Calculate pressure surface once and cache it
            P_mesh = van_der_waals_equation(V_mesh, T_mesh, a_val, b_val, n_val)
            
            # Cache the results
            st.session_state.surface_cache_key = cache_key
            st.session_state.V_range = V_range
            st.session_state.T_range = T_range
            st.session_state.V_mesh = V_mesh
            st.session_state.T_mesh = T_mesh
            st.session_state.P_mesh = P_mesh
    else:
        # Use cached data
        V_range = st.session_state.V_range
        T_range = st.session_state.T_range
        V_mesh = st.session_state.V_mesh
        T_mesh = st.session_state.T_mesh
        P_mesh = st.session_state.P_mesh

    # Calculate pressure slice for current temperature (fast operation)
    V_slice_range = V_range
    P_slice = van_der_waals_equation(V_slice_range, T_slice_val, a_val, b_val, n_val)

    # Create two columns for the plots
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("3D PVT Oberfläche")
        
        # Create 3D plot with cached data
        fig_3d = plt.figure(figsize=(10, 8))
        ax3d = fig_3d.add_subplot(111, projection='3d')
        
        # Create 3D surface plot with reduced complexity for speed
        surface = ax3d.plot_surface(
            V_mesh, T_mesh, P_mesh, cmap='coolwarm', alpha=0.7,
            linewidth=0, antialiased=False, rcount=50, ccount=50  # Reduced surface resolution
        )
        
        # Add temperature slice line - ensure it's always visible on top
        valid_mask = ~np.isnan(P_slice)
        if np.any(valid_mask):
            ax3d.plot(V_slice_range[valid_mask],
                     np.full_like(V_slice_range[valid_mask], T_slice_val),
                     P_slice[valid_mask], 'red', linewidth=3, 
                     zorder=1000, alpha=1.0)
        
        # Set labels and limits
        ax3d.set_xlabel('Volumen (V) [L]')
        ax3d.set_ylabel('Temperatur (T) [K]')
        ax3d.set_zlabel('Druck (P) [Pa]')
        ax3d.set_zlim(0, p_max_val)
        ax3d.set_title(f'Van der Waals: a={a_val:.0f}, b={b_val:.3f}, n={n_val:.1f}')
        
        st.pyplot(fig_3d, clear_figure=True)  # Clear figure for memory optimization
    
    with col2:
        st.subheader("P-V Isotherme (2D Schnitt)")
        
        # Create 2D plot (this is already fast)
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
        
        st.pyplot(fig_2d, clear_figure=True)  # Clear figure for memory optimization
    
    # Display animation if available
    if 'animation_html' in st.session_state:
        st.markdown("---")
        st.subheader("🎬 Temperatur-Animation: 3D + 2D Ansichten")
        st.markdown("**Links: 3D Van der Waals Oberfläche | Rechts: 2D P-V Isotherme**")
        st.markdown("Die rote Linie bewegt sich durch verschiedene Temperaturen (200K bis 400K)")
        
        # Display animation with larger fixed height
        st.components.v1.html(st.session_state.animation_html, height=1000, scrolling=False)
    
    # Add information section
    st.markdown("---")
    st.subheader("Anleitung")
    st.markdown("""
    - **Temperatur-Schnitt**: Verwende den Schieberegler, um die Temperatur der roten Linie im 3D-Plot zu ändern
    - **Parameter a & b**: Gib präzise Van der Waals Parameter für verschiedene Gase ein
    - **Stoffmenge n**: Passe die Anzahl der Mol an
    - **P-Achse Max**: Ändere den maximalen Druck für bessere Visualisierung
    - **🎬 Animation generieren**: Erstellt eine 3D-Animation der sich bewegenden Temperaturlinie
    
    **Animations-Features:**
    - **Doppelte Ansicht**: 3D Van der Waals Oberfläche + 2D P-V Isotherme gleichzeitig
    - **Synchronisiert**: Beide Animationen zeigen die gleiche bewegliche Temperaturlinie  
    - **Temperatur-Sweep**: Rote Linie bewegt sich von 200K bis 400K
    - **Einstellbare Qualität**: 20-80 Frames für verschiedene Glätte-Level
    - **Browser-optimiert**: Stabile HTML5-Animation ohne Server-Belastung
    - **Dauerschleife**: Animation läuft kontinuierlich mit Kontrollen
    
    **Physikalische Bedeutung:**
    - Die **rote Linie** zeigt P-V Verhalten bei konstanter Temperatur (Isotherme)
    - **Niedrige T**: Steile Kurven, starke Druckanstiege
    - **Hohe T**: Flachere Kurven, mehr ideales Gasverhalten
    
    **Aktuelle CO₂ Parameter:**
    - a = 364,0 Pa·L²/mol² (zwischenmolekulare Anziehung)
    - b = 0,04267 L/mol (Molekülvolumen)
    """)
    
    # Close matplotlib figures to prevent memory leaks
    plt.close('all')

if __name__ == "__main__":
    create_3d_pbt_diagram()
