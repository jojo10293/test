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
    R = 8.314 * 10**3  # Universal gas constant [L·kPa/(mol·K)]

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
    
    # Set up figure with fixed, larger size and better spacing
    fig = plt.figure(figsize=(16*0.8, 8*0.8))
    
    # 3D subplot with more space allocation
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
    
    # Adjust layout to prevent cutoff of both 3D plot and titles
    plt.subplots_adjust(left=0.08, right=0.92, top=0.85, bottom=0.15, wspace=0.35)
    
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
    
    # Add logo in upper left corner
    try:
        with open("logo_owks_white.svg", "r", encoding='utf-8') as f:
            logo_svg = f.read()
        
        # Clean the SVG content thoroughly
        import re
        
        # Remove any potential HTML artifacts or comments
        logo_svg = re.sub(r'<!--.*?-->', '', logo_svg, flags=re.DOTALL)
        logo_svg = logo_svg.strip()
        
        # Ensure it's a proper SVG
        if not logo_svg.startswith('<?xml') and not logo_svg.startswith('<svg'):
            raise ValueError("Invalid SVG format")
        
        # Extract just the SVG part if there's XML declaration
        svg_match = re.search(r'<svg.*?</svg>', logo_svg, flags=re.DOTALL)
        if svg_match:
            logo_svg = svg_match.group(0)
        
        # Display logo using st.image with better control
        # First try to display as HTML
        st.markdown(f"""
        <div style="width: 400px; height: auto; margin-bottom: 20px;">
            {logo_svg}
        </div>
        """, unsafe_allow_html=True)
        
        # Title below the logo
        st.title("3D PVT Diagramm - Van der Waals Zustandsgleichung")
            
    except (FileNotFoundError, ValueError, UnicodeDecodeError, Exception) as e:
        # If all else fails, just show the title without logo
        st.title("3D PVT Diagramm - Van der Waals Zustandsgleichung")
        st.info("Logo konnte nicht geladen werden.")
    
    # Sidebar for controls
    st.sidebar.header("Parameter")
    
    # Initial parameters for Van der Waals equation - CO₂ values
    T_init = 300
    a_init = 364.0 * 10**3  # Van der Waals parameter a for CO₂ [Pa·L²/mol²]
    b_init = 0.04267  # Van der Waals parameter b for CO₂ [L/mol]
    n_init = 1.0  # Changed to float
    P_max_init = 13000000  # Reset to reasonable default
    
    # Streamlit controls in sidebar
    T_slice_val = st.sidebar.slider('Temperatur-Schnitt (K)', min_value=200, max_value=400, value=T_init, step=1, key="temp_slider")
    
    # Text inputs for a and b parameters
    col1, col2 = st.sidebar.columns(2)
    with col1:
        a_val = st.number_input('Parameter a', value=a_init, format="%.1f", key="param_a_input")
    with col2:
        b_val = st.number_input('Parameter b', value=b_init, format="%.5f", key="param_b_input")
    
    n_val = st.sidebar.slider('Stoffmenge n (mol)', min_value=0.1, max_value=5.0, value=n_init, step=0.1, key="n_slider")
    p_max_val = st.sidebar.slider('P-Achse Max (Pa)', min_value=1000, max_value=15000000, value=P_max_init, step=1000, key="p_max_slider")
    
    # Animation controls
    st.sidebar.markdown("---")
    st.sidebar.subheader("Animation")
    
    # Animation button and speed control
    col_anim1, col_anim2 = st.sidebar.columns(2)
    with col_anim1:
        generate_animation = st.button("🎬 Animation generieren", key="generate_anim")
    with col_anim2:
        animation_frames = st.selectbox("Frames", [20, 40, 60, 80, 120], index=2, key="frames")

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
            try:
                # Create volume and temperature ranges with optimized resolution
                V_range = np.linspace(0.00001, 1, 150)  # Reduced from 200 for speed
                T_range = np.linspace(200, 400, 150)
                V_mesh, T_mesh = np.meshgrid(V_range, T_range)

                # Calculate pressure surface once and cache it
                P_mesh = van_der_waals_equation(V_mesh, T_mesh, a_val, b_val, n_val)
                
                # Handle potential numerical issues
                #P_mesh = np.where(np.isfinite(P_mesh), P_mesh, np.nan)
                #P_mesh = np.where(P_mesh > 0, P_mesh, np.nan)  # Mask non-physical pressures
                
                # Cache the results
                st.session_state.surface_cache_key = cache_key
                st.session_state.V_range = V_range
                st.session_state.T_range = T_range
                st.session_state.V_mesh = V_mesh
                st.session_state.T_mesh = T_mesh
                st.session_state.P_mesh = P_mesh
                
            except Exception as e:
                st.error(f"Fehler bei der Berechnung: {str(e)}")
                st.error("Bitte Parameter anpassen!")
                return  # Exit function if calculation fails
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
    
    # Handle numerical issues in slice calculation
    P_slice = np.where(np.isfinite(P_slice), P_slice, np.nan)
    #P_slice = np.where(P_slice > 0, P_slice, np.nan)

    # Create combined figure with both plots (using same design as animation)
    st.subheader("Van der Waals PVT Diagramm: 3D Oberfläche & 2D Isotherme")
    
    # Create single figure with smaller size for more white space
    fig_combined = plt.figure(figsize=(14, 7))  # Reduced from (16*0.8, 8*0.8) = (12.8, 6.4)
    
    # 3D subplot (left side) - same settings as animation
    ax3d = fig_combined.add_subplot(121, projection='3d')
    
    # Create 3D surface plot with same settings as animation
    surface = ax3d.plot_surface(
        V_mesh, T_mesh, P_mesh, cmap='coolwarm', alpha=0.6,  # Same alpha as animation
        linewidth=0, antialiased=False  # Same settings as animation
    )
    
    # Add temperature slice line with same settings as animation
    valid_mask = ~np.isnan(P_slice)
    if np.any(valid_mask):
        ax3d.plot(V_slice_range[valid_mask],
                 np.full_like(V_slice_range[valid_mask], T_slice_val),
                 P_slice[valid_mask], 'red', linewidth=3, 
                 zorder=1000, alpha=1.0)  # Same settings as animation
    
    # Set 3D plot labels and limits - same as animation
    ax3d.set_xlabel('Volumen (V) [L]', fontsize=10)  # Same fontsize as animation
    ax3d.set_ylabel('Temperatur (T) [K]', fontsize=10)  # Same fontsize as animation
    ax3d.set_zlabel('Druck (P) [Pa]', fontsize=10)  # Same fontsize as animation
    ax3d.set_xlim(0, 1)  # Same limits as animation
    ax3d.set_ylim(200, 400)  # Same limits as animation
    ax3d.set_zlim(0, p_max_val)
    ax3d.set_title(f'Van der Waals 3D\nT = {T_slice_val:.0f} K', fontsize=12)  # Same title style as animation
    
    # 2D subplot (right side) - same settings as animation
    ax2d = fig_combined.add_subplot(122)
    
    # Plot 2D slice with same settings as animation
    valid_mask_2d = ~np.isnan(P_slice)
    if np.any(valid_mask_2d):
        ax2d.plot(V_slice_range[valid_mask_2d], P_slice[valid_mask_2d], 'red', linewidth=3)  # Same as animation
    
    # Set 2D plot labels and limits - same as animation
    ax2d.set_xlabel('Volumen (V) [L]', fontsize=12)  # Same fontsize as animation
    ax2d.set_ylabel('Druck (P) [Pa]', fontsize=12)  # Same fontsize as animation
    ax2d.set_xlim(0, 1)  # Same limits as animation
    ax2d.set_ylim(0, p_max_val)
    ax2d.set_title(f'P-V Isotherme\nT = {T_slice_val:.0f} K', fontsize=12)  # Same title style as animation
    ax2d.grid(True, alpha=0.3)  # Same grid settings as animation
    
    # Use more generous layout with bigger padding for less cramped appearance
    plt.subplots_adjust(left=0.12, right=0.88, top=0.82, bottom=0.18, wspace=0.4)
    
    # Display the combined figure
    st.pyplot(fig_combined, clear_figure=True)
    
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
    Mit diesem Tool kannst du das Verhalten eines realen Gases nach der Van der Waals Gleichung interaktiv erkunden. 
    Über den Temperatur-Schnitt-Schieberegler lässt sich die rote Isotherme im 3D-Plot direkt verändern. 
    Für einen anschaulichen Überblick kannst du ausserdem eine Animation generieren, in der die rote Linie flüssig durch verschiedene Temperaturen (200K bis 400K) läuft – so siehst du, wie sich die P-V-Isothermen mit steigender Temperatur verändern. 
    Passe die Van der Waals Parameter (a, b), die Stoffmenge n und den maximalen Druck nach Bedarf an, um unterschiedliche Gase oder Bedingungen zu simulieren.
    Die Animation zeigt dabei synchron die Entwicklung im 3D-Diagramm und als 2D-Schnitt.
    """)
    
    # Close matplotlib figures to prevent memory leaks
    plt.close('all')

if __name__ == "__main__":
    create_3d_pbt_diagram()
