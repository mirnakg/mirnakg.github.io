import streamlit as st
import os
import tempfile
import subprocess
import urllib.request
import stat
from math import sqrt

st.set_page_config(
    page_title="Vina Plug & Play",
    page_icon="🧬",
    layout="wide",
)

# --- Custom CSS to match portfolio aesthetic ---
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;700&display=swap');

html, body, [class*="css"] {
    font-family: 'IBM Plex Mono', monospace;
}
.stApp {
    background-color: #f0ece3;
}
h1, h2, h3 {
    font-family: 'IBM Plex Mono', monospace;
    letter-spacing: 0.05em;
}
.stButton > button {
    background-color: #2a52cc;
    color: white;
    border: none;
    font-family: 'IBM Plex Mono', monospace;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    font-weight: 700;
    padding: 0.6em 2em;
    transition: opacity 0.2s;
}
.stButton > button:hover {
    background-color: #1e3fa0;
    color: white;
}
div[data-testid="stSidebar"] {
    background-color: #e8e4db;
    border-right: 1.5px dashed #c0bbb0;
}
.step-header {
    font-size: 0.7em;
    letter-spacing: 0.25em;
    text-transform: uppercase;
    color: #2a52cc;
    margin-bottom: 0.5rem;
    font-weight: 700;
}
</style>
""", unsafe_allow_html=True)


# --- Vina binary setup ---
@st.cache_resource
def get_vina_binary():
    """Download and cache the AutoDock Vina binary."""
    vina_dir = os.path.join(tempfile.gettempdir(), "vina_bin")
    os.makedirs(vina_dir, exist_ok=True)
    vina_path = os.path.join(vina_dir, "vina_1.2.5_linux_x86_64")

    if not os.path.exists(vina_path):
        url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
        urllib.request.urlretrieve(url, vina_path)
        os.chmod(vina_path, os.stat(vina_path).st_mode | stat.S_IEXEC)

    return vina_path


# --- Helper Functions (all conversions via Open Babel CLI) ---

def smiles_to_pdbqt(smiles, output_path):
    """Convert SMILES to 3D PDBQT using Open Babel."""
    # Step 1: SMILES -> 3D SDF (with hydrogen addition and geometry optimization)
    sdf_path = output_path.replace(".pdbqt", ".sdf")
    result = subprocess.run(
        ["obabel", f"-:{smiles}", "-osdf", "-O", sdf_path,
         "--gen3d", "-h"],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not os.path.exists(sdf_path):
        raise ValueError(f"Could not generate 3D structure from SMILES: {result.stderr}")

    # Step 2: SDF -> PDBQT
    result = subprocess.run(
        ["obabel", sdf_path, "-O", output_path, "-xh"],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not os.path.exists(output_path):
        raise ValueError(f"Could not convert to PDBQT: {result.stderr}")

    return output_path


def fetch_pdb(pdb_id, output_path):
    """Download PDB file from RCSB."""
    import requests
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
    r = requests.get(url)
    if r.status_code != 200:
        raise ValueError(f"Could not download PDB ID '{pdb_id}'. Check the ID and try again.")
    with open(output_path, "w") as f:
        f.write(r.text)
    return output_path


def detect_ligands(pdb_path):
    """Detect hetero residues (non-water) in PDB file."""
    ligands = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                resn = line[17:20].strip()
                if resn != "HOH":
                    ligands.add(resn)
    return ligands


def get_ligand_center(pdb_path, ligand_resn):
    """Calculate center of mass of a ligand residue."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                resn = line[17:20].strip()
                if resn == ligand_resn:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
    if not coords:
        raise ValueError(f"Could not find ligand '{ligand_resn}' in the PDB file.")
    center = (
        sum(c[0] for c in coords) / len(coords),
        sum(c[1] for c in coords) / len(coords),
        sum(c[2] for c in coords) / len(coords),
    )
    return center


def clean_receptor(pdb_path, output_path, ligands_to_remove, keep_metals=True,
                   mode="full", pocket_center=None, pocket_radius=15):
    """Clean receptor PDB: remove ligands, optionally trim to pocket."""
    metal_ions = {"ZN", "MG", "FE", "CA", "MN", "CU", "CO", "NI"}

    def dist(a, b):
        return sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

    with open(pdb_path) as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                if mode == "full":
                    outfile.write(line)
                elif mode == "pocket" and pocket_center:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if dist((x, y, z), pocket_center) <= pocket_radius:
                        outfile.write(line)
            elif line.startswith("HETATM"):
                resn = line[17:20].strip()
                if resn in ligands_to_remove:
                    continue
                if keep_metals and resn in metal_ions:
                    outfile.write(line)

    return output_path


def convert_to_pdbqt(pdb_path, pdbqt_path):
    """Convert PDB to PDBQT using Open Babel."""
    result = subprocess.run(
        ["obabel", pdb_path, "-O", pdbqt_path, "-xr", "-xp", "-h"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel conversion failed: {result.stderr}")
    return pdbqt_path


def run_docking(vina_bin, receptor_pdbqt, ligand_pdbqt, center, box_size,
                output_path, exhaustiveness=32, n_poses=10):
    """Run AutoDock Vina docking via CLI binary."""
    cmd = [
        vina_bin,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--out", output_path,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"Vina docking failed:\n{result.stderr}")

    # Parse energies from the output PDBQT file
    energies = []
    with open(output_path) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                energy = float(parts[3])
                rmsd_lb = float(parts[4])
                rmsd_ub = float(parts[5])
                energies.append((energy, rmsd_lb, rmsd_ub))

    return energies, output_path


def extract_poses_to_pdb(poses_pdbqt_path, output_dir, best_only=False):
    """Extract docked poses from PDBQT to individual PDB files using Open Babel."""
    # Split multi-model PDBQT into individual models
    models = []
    current_model = []
    with open(poses_pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = []
            elif line.startswith("ENDMDL"):
                models.append("".join(current_model))
            else:
                current_model.append(line)

    if not models and current_model:
        models.append("".join(current_model))

    pose_files = []
    for i, model_content in enumerate(models):
        single_pdbqt = os.path.join(output_dir, f"pose_{i+1}.pdbqt")
        pose_pdb = os.path.join(output_dir, f"pose_{i+1}.pdb")

        with open(single_pdbqt, "w") as f:
            f.write(model_content)

        subprocess.run(
            ["obabel", single_pdbqt, "-O", pose_pdb],
            capture_output=True, text=True
        )

        if os.path.exists(pose_pdb):
            pose_files.append(pose_pdb)

        if best_only:
            break

    return pose_files


# ============================================================
# APP LAYOUT
# ============================================================

st.markdown("# Vina Plug & Play")
st.markdown("**Interactive molecular docking with AutoDock Vina** — no programming required.")
st.markdown("---")

# Sidebar
with st.sidebar:
    st.markdown("### How it works")
    st.markdown("""
1. **Provide a ligand** via SMILES string or .pdbqt file upload
2. **Provide a receptor** via PDB ID or .pdb file upload
3. **Configure** the binding pocket and cleaning options
4. **Run docking** and view ranked binding poses
    """)
    st.markdown("---")
    st.markdown(
        '<p style="font-size:0.7em; letter-spacing:0.15em; text-transform:uppercase;">'
        'Built by <a href="https://github.com/mirnakg" style="color:#2a52cc;">Mirna Kheir Gouda</a>'
        '</p>',
        unsafe_allow_html=True
    )

# Temp directory for session
if "work_dir" not in st.session_state:
    st.session_state.work_dir = tempfile.mkdtemp()
work_dir = st.session_state.work_dir

# --- STEP 1: LIGAND ---
st.markdown('<div class="step-header">Step 1 — Ligand</div>', unsafe_allow_html=True)
ligand_method = st.radio("How would you like to provide your ligand?",
                         ["SMILES string", "Upload .pdbqt file"], horizontal=True)

ligand_ready = False
ligand_pdbqt = None

if ligand_method == "SMILES string":
    smiles = st.text_input("Enter SMILES string", placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O")
    if smiles:
        ligand_ready = True
else:
    ligand_upload = st.file_uploader("Upload ligand .pdbqt", type=["pdbqt"])
    if ligand_upload:
        ligand_path = os.path.join(work_dir, ligand_upload.name)
        with open(ligand_path, "wb") as f:
            f.write(ligand_upload.read())
        ligand_pdbqt = ligand_path
        ligand_ready = True

st.markdown("---")

# --- STEP 2: RECEPTOR ---
st.markdown('<div class="step-header">Step 2 — Receptor</div>', unsafe_allow_html=True)
receptor_method = st.radio("How would you like to provide your receptor?",
                           ["PDB ID (fetch from RCSB)", "Upload .pdb file"], horizontal=True)

receptor_ready = False
pdb_path = None

if receptor_method == "PDB ID (fetch from RCSB)":
    pdb_id = st.text_input("Enter PDB ID", placeholder="e.g. 4UOX")
    if pdb_id:
        receptor_ready = True
else:
    receptor_upload = st.file_uploader("Upload receptor .pdb", type=["pdb"])
    if receptor_upload:
        pdb_path = os.path.join(work_dir, receptor_upload.name)
        with open(pdb_path, "wb") as f:
            f.write(receptor_upload.read())
        receptor_ready = True

st.markdown("---")

# --- STEP 3: CONFIGURATION ---
st.markdown('<div class="step-header">Step 3 — Configuration</div>', unsafe_allow_html=True)

col1, col2 = st.columns(2)

with col1:
    keep_metals = st.checkbox("Keep metal ions", value=True)
    receptor_mode = st.selectbox("Receptor mode", ["Full", "Pocket"])

with col2:
    pocket_method = st.selectbox("Pocket center definition",
                                  ["From ligand residue", "Manual coordinates"])
    pocket_radius = st.number_input("Pocket radius (Angstrom)", value=15.0, min_value=1.0, max_value=50.0)

if pocket_method == "From ligand residue":
    pocket_resn = st.text_input("Ligand residue name for pocket center", placeholder="e.g. PLP").strip().upper()
else:
    mcol1, mcol2, mcol3 = st.columns(3)
    manual_x = mcol1.number_input("Center X", value=0.0, format="%.2f")
    manual_y = mcol2.number_input("Center Y", value=0.0, format="%.2f")
    manual_z = mcol3.number_input("Center Z", value=0.0, format="%.2f")

col_ex, col_np = st.columns(2)
with col_ex:
    exhaustiveness = st.slider("Exhaustiveness", min_value=8, max_value=64, value=32, step=8)
with col_np:
    n_poses = st.slider("Number of poses", min_value=1, max_value=20, value=10)

st.markdown("---")

# --- STEP 4: RUN ---
st.markdown('<div class="step-header">Step 4 — Run Docking</div>', unsafe_allow_html=True)

can_run = ligand_ready and receptor_ready
if pocket_method == "From ligand residue" and not pocket_resn:
    can_run = False

if st.button("Run Docking", disabled=not can_run, type="primary"):
    progress = st.progress(0, text="Initializing...")

    try:
        # 0. Get vina binary
        progress.progress(5, text="Setting up AutoDock Vina...")
        vina_bin = get_vina_binary()

        # 1. Prepare ligand
        progress.progress(10, text="Preparing ligand...")
        if ligand_pdbqt is None:
            ligand_pdbqt = os.path.join(work_dir, "ligand.pdbqt")
            smiles_to_pdbqt(smiles, ligand_pdbqt)
            st.success("Ligand generated from SMILES")

        # 2. Prepare receptor
        progress.progress(25, text="Preparing receptor...")
        if pdb_path is None:
            pdb_path = os.path.join(work_dir, f"{pdb_id.lower()}.pdb")
            fetch_pdb(pdb_id, pdb_path)
            st.success(f"Downloaded PDB: {pdb_id}")

        # Detect and display ligands
        detected_ligands = detect_ligands(pdb_path)
        if detected_ligands:
            st.info(f"Detected hetero residues: {', '.join(sorted(detected_ligands))}")

        # 3. Determine pocket center
        progress.progress(35, text="Defining binding pocket...")
        if pocket_method == "From ligand residue":
            center = list(get_ligand_center(pdb_path, pocket_resn))
            st.success(f"Pocket center from {pocket_resn}: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
        else:
            center = [manual_x, manual_y, manual_z]
            st.success(f"Manual pocket center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")

        # 4. Clean receptor
        progress.progress(45, text="Cleaning receptor...")
        clean_pdb = os.path.join(work_dir, "receptor_clean.pdb")
        mode = "pocket" if receptor_mode == "Pocket" else "full"
        clean_receptor(pdb_path, clean_pdb, detected_ligands, keep_metals=keep_metals,
                       mode=mode, pocket_center=tuple(center), pocket_radius=pocket_radius)

        # 5. Convert to PDBQT
        progress.progress(55, text="Converting receptor to PDBQT...")
        receptor_pdbqt = os.path.join(work_dir, "receptor_clean.pdbqt")
        convert_to_pdbqt(clean_pdb, receptor_pdbqt)

        # 6. Run docking
        box_dim = 2 * pocket_radius
        box_size = [box_dim, box_dim, box_dim]
        progress.progress(65, text="Running AutoDock Vina (this may take a few minutes)...")
        poses_path = os.path.join(work_dir, "docked_poses.pdbqt")
        energies, poses_path = run_docking(vina_bin, receptor_pdbqt, ligand_pdbqt,
                                           center, box_size, poses_path,
                                           exhaustiveness=exhaustiveness, n_poses=n_poses)

        progress.progress(90, text="Extracting poses...")

        # 7. Extract poses to PDB
        pose_files = extract_poses_to_pdb(poses_path, work_dir, best_only=False)

        progress.progress(100, text="Done!")

        # --- RESULTS ---
        st.markdown("---")
        st.markdown('<div class="step-header">Results</div>', unsafe_allow_html=True)

        # Energy table
        st.markdown("**Binding Energies (kcal/mol)**")
        import pandas as pd
        energy_data = []
        for i, e in enumerate(energies):
            energy_data.append({
                "Pose": i + 1,
                "Affinity": f"{e[0]:.2f}",
                "RMSD l.b.": f"{e[1]:.2f}",
                "RMSD u.b.": f"{e[2]:.2f}",
            })
        df = pd.DataFrame(energy_data)
        st.dataframe(df, use_container_width=True, hide_index=True)

        st.success(f"Best binding affinity: **{energies[0][0]:.2f} kcal/mol**")

        # Download buttons
        st.markdown("**Download Results**")
        dcol1, dcol2 = st.columns(2)

        with dcol1:
            with open(poses_path, "r") as f:
                st.download_button("Download all poses (.pdbqt)", f.read(),
                                   file_name="docked_poses.pdbqt", mime="chemical/x-pdbqt")
        with dcol2:
            if pose_files:
                with open(pose_files[0], "r") as f:
                    st.download_button("Download best pose (.pdb)", f.read(),
                                       file_name="best_pose.pdb", mime="chemical/x-pdb")

    except Exception as e:
        progress.empty()
        st.error(f"Error: {str(e)}")

elif not can_run:
    missing = []
    if not ligand_ready:
        missing.append("ligand")
    if not receptor_ready:
        missing.append("receptor")
    if pocket_method == "From ligand residue" and not pocket_resn:
        missing.append("pocket residue name")
    st.caption(f"Provide {', '.join(missing)} to enable docking.")
