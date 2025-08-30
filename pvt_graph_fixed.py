import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import streamlit as st
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
    
    # Animation controls
    st.sidebar.markdown("---")
    st.sidebar.subheader("Animation")
    
    # Animation button and speed control
    col_anim1, col_anim2 = st.sidebar.columns(2)
    with col_anim1:
        animate_button = st.button("🎬 Animation starten", key="animate")
    with col_anim2:
        animation_speed = st.selectbox("Geschwindigkeit", [0.1, 0.2, 0.5, 1.0], index=2, key="speed")
    
    # Initialize session state for animation
    if 'animating' not in st.session_state:
        st.session_state.animating = False
    if 'current_temp' not in st.session_state:
        st.session_state.current_temp = T_init
    if 'temp_direction' not in st.session_state:
        st.session_state.temp_direction = 1  # 1 for increasing, -1 for decreasing
    
    # Handle animation
    if animate_button:
        st.session_state.animating = not st.session_state.animating
    
    # Animation logic
    if st.session_state.animating:
        # Update temperature for animation
        temp_step = 5 * st.session_state.temp_direction
        new_temp = st.session_state.current_temp + temp_step
        
        # Reverse direction at boundaries
        if new_temp >= 400:
            st.session_state.temp_direction = -1
            new_temp = 400
        elif new_temp <= 200:
            st.session_state.temp_direction = 1
            new_temp = 200
            
        st.session_state.current_temp = new_temp
        T_slice_val = new_temp
        
        # Add animation status indicator
        st.sidebar.success(f"🎬 Animiere bei T = {T_slice_val:.0f} K")
        st.sidebar.info("Klicke 'Animation starten' erneut zum Stoppen")
        
        # Auto-refresh for animation
        time.sleep(animation_speed)
        st.rerun()
    else:
        # Use slider value when not animating
        st.session_state.current_temp = T_slice_val
    
    # Text inputs for a and b parameters
    col1, col2 = st.sidebar.columns(2)
    with col1:
        a_val = st.number_input('Parameter a', value=a_init, format="%.1f", key="param_a")
    with col2:
        b_val = st.number_input('Parameter b', value=b_init, format="%.5f", key="param_b")
    
    n_val = st.sidebar.slider('Stoffmenge n (mol)', min_value=0.1, max_value=5.0, value=n_init, step=0.1)
    p_max_val = st.sidebar.slider('P-Achse Max (Pa)', min_value=1000, max_value=100000, value=P_max_init, step=1000)
    
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
        animation_status = " (ANIMIEREND)" if st.session_state.animating else ""
        ax3d.set_title(f'Van der Waals: a={a_val:.0f}, b={b_val:.3f}, n={n_val:.1f}{animation_status}')
        
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
        
        # Add animation indicator to 2D title as well
        animation_status = " (ANIMIEREND)" if st.session_state.animating else ""
        ax2d.set_title(f'P-V Isotherme bei T = {T_slice_val:.0f} K{animation_status}')
        ax2d.grid(True, alpha=0.3)
        
        st.pyplot(fig_2d)
    
    # Add information section
    st.markdown("---")
    st.subheader("Anleitung")
    st.markdown("""
    - **Temperatur-Schnitt**: Verwende den Schieberegler, um die Temperatur der roten Linie im 3D-Plot zu ändern
    - **Parameter a & b**: Gib präzise Van der Waals Parameter für verschiedene Gase ein
    - **Stoffmenge n**: Passe die Anzahl der Mol an
    - **P-Achse Max**: Ändere den maximalen Druck für bessere Visualisierung
    - **🎬 Animation**: Klicke "Animation starten" um automatisch durch die Temperaturen zu fahren
    
    **Animations-Features:**
    - Läuft automatisch von 200K bis 400K und zurück
    - Einstellbare Geschwindigkeit (0,1s bis 1,0s pro Frame)
    - Erneut klicken zum Stoppen
    - Zeigt dynamisches Van der Waals Verhalten über den Temperaturbereich
    
    **Aktuelle CO₂ Parameter:**
    - a = 364,0 Pa·L²/mol² (zwischenmolekulare Anziehung)
    - b = 0,04267 L/mol (Molekülvolumen)
    """)
    
    # Close matplotlib figures to prevent memory leaks
    plt.close('all')

if __name__ == "__main__":
    create_3d_pbt_diagram()
