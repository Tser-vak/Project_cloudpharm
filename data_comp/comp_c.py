import pandas as pd
import os
import re

class DockingComparer:
    def __init__(self, bpp_path, old_path, output_dir, boltz_path=None):
        self.bpp_path = bpp_path
        self.old_path = old_path
        self.output_dir = output_dir
        self.boltz_path = boltz_path
        
        self.df_bpp = None
        self.df_old = None
        self.df_merged = None
        
        self.score_cols = ['Vina_Affinity(kcal/mol)', 'CNN_Pose_Score']

    def load_data(self):
        self.df_bpp = pd.read_csv(self.bpp_path)
        self.df_old = pd.read_csv(self.old_path)
        
        # Handle potential "FAILED" strings by forcing columns to numeric, coercing errors to NaN
        for col in self.score_cols:
            self.df_bpp[col] = pd.to_numeric(self.df_bpp[col], errors='coerce')
            self.df_old[col] = pd.to_numeric(self.df_old[col], errors='coerce')
            
        # Merge the dataframes on Ligand_Name
        self.df_merged = pd.merge(
            self.df_old[['Ligand_Name'] + self.score_cols], 
            self.df_bpp[['Ligand_Name'] + self.score_cols],
            on='Ligand_Name', 
            suffixes=('_old', '_new')
        )
        
        # Calculate differences -> Vina difference of absolute values (magnitude)
        self.df_merged['Vina_diff_new_minus_old'] = self.df_merged['Vina_Affinity(kcal/mol)_new'].abs() - self.df_merged['Vina_Affinity(kcal/mol)_old'].abs()
        self.df_merged['CNN_diff_new_minus_old'] = self.df_merged['CNN_Pose_Score_new'] - self.df_merged['CNN_Pose_Score_old']
        
        # Calculate winners
        self.df_merged['Better_Vina'] = self.df_merged.apply(self._get_better_vina, axis=1)
        self.df_merged['Better_CNN'] = self.df_merged.apply(self._get_better_cnn, axis=1)
        self.df_merged['Result_Match'] = self.df_merged.apply(self._compare_overall_match, axis=1)

    def _get_better_vina(self, row):
        if pd.isna(row['Vina_diff_new_minus_old']):
            return "Unknown"
        if row['Vina_Affinity(kcal/mol)_new'] < row['Vina_Affinity(kcal/mol)_old']:
            return "New (BPP)"
        elif row['Vina_Affinity(kcal/mol)_new'] > row['Vina_Affinity(kcal/mol)_old']:
            return "Old"
        return "Same"

    def _get_better_cnn(self, row):
        if pd.isna(row['CNN_diff_new_minus_old']):
            return "Unknown"
        if row['CNN_Pose_Score_new'] > row['CNN_Pose_Score_old']:
            return "New (BPP)"
        elif row['CNN_Pose_Score_new'] < row['CNN_Pose_Score_old']:
            return "Old"
        return "Same"

    def _compare_overall_match(self, row):
        if pd.isna(row['Vina_diff_new_minus_old']) or pd.isna(row['CNN_diff_new_minus_old']):
            return "Unknown/Failed"
        if row['Vina_diff_new_minus_old'] == 0 and row['CNN_diff_new_minus_old'] == 0:
            return "Exactly Same"
        return "Different"

    def generate_overall_summary(self):
        output_path = os.path.join(self.output_dir, 'comparison_results.csv')
        self.df_merged.to_csv(output_path, index=False)
        
        total_complexes = len(self.df_merged)
        exactly_same = len(self.df_merged[self.df_merged['Result_Match'] == 'Exactly Same'])
        different = len(self.df_merged[self.df_merged['Result_Match'] == 'Different'])
        unknown_failed = len(self.df_merged[self.df_merged['Result_Match'] == 'Unknown/Failed'])
        
        vina_better_old = len(self.df_merged[self.df_merged['Better_Vina'] == 'Old'])
        vina_better_new = len(self.df_merged[self.df_merged['Better_Vina'] == 'New (BPP)'])
        vina_same = len(self.df_merged[self.df_merged['Better_Vina'] == 'Same'])
        
        cnn_better_old = len(self.df_merged[self.df_merged['Better_CNN'] == 'Old'])
        cnn_better_new = len(self.df_merged[self.df_merged['Better_CNN'] == 'New (BPP)'])
        cnn_same = len(self.df_merged[self.df_merged['Better_CNN'] == 'Same'])
        
        print("--- DOCKING RESULTS COMPARISON SUMMARY ---")
        print(f"Total complexes compared: {total_complexes}\n")
        print("Match Results:")
        print(f"  Complexes with exactly the same results: {exactly_same}")
        print(f"  Complexes with different results: {different}")
        print(f"  Complexes with errors/failed docking: {unknown_failed}\n")
        print("Vina Affinity (Lower is Better):")
        print(f"  Better in OLD: {vina_better_old}")
        print(f"  Better in NEW (BPP): {vina_better_new}")
        print(f"  Same in both: {vina_same}\n")
        print("CNN Pose Score (Higher is Better):")
        print(f"  Better in OLD: {cnn_better_old}")
        print(f"  Better in NEW (BPP): {cnn_better_new}")
        print(f"  Same in both: {cnn_same}\n")
        print(f"Detailed comparison saved to: '{output_path}'")

    def _get_top_20_percent(self, df, sort_col, ascending=True):
        top_20_count = int(len(df) * 0.20)
        return df.sort_values(sort_col, ascending=ascending).head(top_20_count)

    def _write_top_results_to_txt(self, df_affinity, df_cnn, version_name, output_path, affinity_col, cnn_col):
        count = len(df_affinity)
        with open(output_path, 'w') as f:
            f.write(f"--- Top 20% ({count} complexes) Analysis for {version_name} ---\n\n")
            
            f.write("--> Top 20% by Vina Affinity (Lower is better):\n")
            f.write(f"    {'Ligand_Name':<35} | {'Affinity':<10} | {'CNN Pose':<10}\n")
            f.write("    " + "-"*60 + "\n")
            for _, row in df_affinity.iterrows():
                f.write(f"    {row['Ligand_Name']:<35} | {row[affinity_col]:<10.2f} | {row[cnn_col]:<10.4f}\n")
            
            f.write("\n--> Top 20% by CNN Pose Score (Higher is better):\n")
            f.write(f"    {'Ligand_Name':<35} | {'CNN Pose':<10} | {'Affinity':<10}\n")
            f.write("    " + "-"*60 + "\n")
            for _, row in df_cnn.iterrows():
                f.write(f"    {row['Ligand_Name']:<35} | {row[cnn_col]:<10.4f} | {row[affinity_col]:<10.2f}\n")
                
        print(f"Top 20% analysis for {version_name} saved to text file: '{output_path}'")

    def _compare_top_20_sets(self, df_old_top, df_new_top, category_name, output_path, winner_col):
        common_ligands = set(df_old_top['Ligand_Name']).intersection(set(df_new_top['Ligand_Name']))
        df_common = self.df_merged[self.df_merged['Ligand_Name'].isin(common_ligands)]
        
        with open(output_path, 'a') as f:
            f.write(f"\n================ CROSS-VERSION COMPARISON: {category_name} ================\n")
            f.write(f"Common complexes found in both Old and New top 20%: {len(common_ligands)}\n")
            
            if len(df_common) > 0:
                f.write(f"\n    {'Ligand_Name':<35} | {'Vina Abs Diff':<13} | {'CNN Diff':<12} | {'Best Score In':<15}\n")
                f.write("    " + "-"*85 + "\n")
                for _, row in df_common.iterrows():
                    f.write(f"    {row['Ligand_Name']:<35} | {row['Vina_diff_new_minus_old']:<13.2f} | {row['CNN_diff_new_minus_old']:<12.4f} | {row[winner_col]:<15}\n")
            else:
                f.write("\n    No common complexes found.\n")
        
        print(f"Comparison for {category_name} results appended to: '{output_path}'")

    def run_top_20_analysis(self):
        # Top sets for Old Dataset
        old_top_affinity = self._get_top_20_percent(self.df_merged, 'Vina_Affinity(kcal/mol)_old', ascending=True)
        old_top_cnn = self._get_top_20_percent(self.df_merged, 'CNN_Pose_Score_old', ascending=False)
        old_txt_path = os.path.join(self.output_dir, 'top_20_old.txt')
        self._write_top_results_to_txt(
            old_top_affinity, old_top_cnn, "OLD DATASET", old_txt_path, 
            'Vina_Affinity(kcal/mol)_old', 'CNN_Pose_Score_old'
        )

        # Top sets for New (BPP) Dataset
        new_top_affinity = self._get_top_20_percent(self.df_merged, 'Vina_Affinity(kcal/mol)_new', ascending=True)
        new_top_cnn = self._get_top_20_percent(self.df_merged, 'CNN_Pose_Score_new', ascending=False)
        new_txt_path = os.path.join(self.output_dir, 'top_20_new.txt')
        self._write_top_results_to_txt(
            new_top_affinity, new_top_cnn, "NEW (BPP) DATASET", new_txt_path, 
            'Vina_Affinity(kcal/mol)_new', 'CNN_Pose_Score_new'
        )

        # Cross-comparison over the common items
        comparison_txt_path = os.path.join(self.output_dir, 'top_20_comparison.txt')
        with open(comparison_txt_path, 'w') as f:
            f.write("--- TOP 20% CROSS-VERSION OVERLAP ANALYSIS ---\n")
        
        self._compare_top_20_sets(
            old_top_affinity, new_top_affinity, 
            "Vina Affinity Top 20%", comparison_txt_path, 'Better_Vina'
        )
        
        self._compare_top_20_sets(
            old_top_cnn, new_top_cnn, 
            "CNN Pose Score Top 20%", comparison_txt_path, 'Better_CNN'
        )

    # ------------------------------------------------------------------
    # Boltz enrichment helpers
    # ------------------------------------------------------------------
    def _parse_top20_txt(self, filename, col_suffix):
        """Generic parser for top_20_new.txt / top_20_old.txt.

        Returns (df_affinity, df_cnn) with columns:
          Ligand_Name,
          Vina_Affinity(kcal/mol){col_suffix},
          CNN_Pose_Score{col_suffix}
        """
        txt_path = os.path.join(self.output_dir, filename)
        if not os.path.exists(txt_path):
            raise FileNotFoundError(
                f"'{filename}' not found at '{txt_path}'.\n"
                "Run run_top_20_analysis() first."
            )

        affinity_col = f'Vina_Affinity(kcal/mol){col_suffix}'
        cnn_col      = f'CNN_Pose_Score{col_suffix}'

        affinity_rows = []
        cnn_rows      = []
        current       = None  # 'affinity' | 'cnn' | None

        row_re = re.compile(
            r'^\s+(\S+)\s*\|\s*(-?[\d.]+)\s*\|\s*(-?[\d.]+)'
        )

        with open(txt_path, 'r', encoding='utf-8') as fh:
            for line in fh:
                stripped = line.strip()

                if 'Top 20% by Vina Affinity' in stripped:
                    current = 'affinity'
                    continue
                if 'Top 20% by CNN Pose Score' in stripped:
                    current = 'cnn'
                    continue
                if not stripped or stripped.startswith('---') or stripped.startswith('Ligand'):
                    continue
                if set(stripped).issubset({'-', ' '}):
                    continue

                m = row_re.match(line)
                if m and current == 'affinity':
                    affinity_rows.append({
                        'Ligand_Name': m.group(1),
                        affinity_col:  float(m.group(2)),
                        cnn_col:       float(m.group(3)),
                    })
                elif m and current == 'cnn':
                    cnn_rows.append({
                        'Ligand_Name': m.group(1),
                        cnn_col:       float(m.group(2)),
                        affinity_col:  float(m.group(3)),
                    })

        return pd.DataFrame(affinity_rows), pd.DataFrame(cnn_rows)

    def _run_boltz_enrichment(self, txt_filename, col_suffix, dataset_label, out_filename):
        """Shared engine for boltz enrichment.
        Reads <txt_filename>, cross-references boltz_dt, writes <out_filename>.
        col_suffix is '_new' or '_old'.
        """
        if self.boltz_path is None:
            print("[WARN] boltz_path not set – skipping Boltz enrichment.")
            return
        if not os.path.exists(self.boltz_path):
            raise FileNotFoundError(f"boltz_dt.csv not found at '{self.boltz_path}'.")

        boltz_cols = ['ligand_id', 'no_2og_binary', 'no_2og_pred_value', 'no_2og_ligand_iptm']
        df_boltz = pd.read_csv(self.boltz_path, usecols=boltz_cols)

        affinity_col = f'Vina_Affinity(kcal/mol){col_suffix}'
        cnn_col      = f'CNN_Pose_Score{col_suffix}'

        df_affinity, df_cnn = self._parse_top20_txt(txt_filename, col_suffix)

        def enrich(df_section, section_label):
            if df_section.empty:
                return pd.DataFrame()
            df = df_section.copy()
            df['_lookup_id'] = df['Ligand_Name'].str.replace(
                r'_model_0_ligand$', '', regex=True
            )
            merged = df.merge(
                df_boltz,
                left_on='_lookup_id',
                right_on='ligand_id',
                how='inner'
            )
            merged = merged.drop(columns=['_lookup_id', 'ligand_id'])
            print(f"  [{section_label}] {len(df)} entries in {txt_filename}  "
                  f"->  {len(merged)} matched in boltz_dt")
            return merged

        enriched_affinity = enrich(df_affinity, 'Vina Affinity')
        enriched_cnn      = enrich(df_cnn,      'CNN Pose Score')

        out_path = os.path.join(self.output_dir, out_filename)
        with open(out_path, 'w', encoding='utf-8') as f:

            def write_section(df_sec, title):
                f.write(f"\n{'='*80}\n")
                f.write(f"  {title}\n")
                f.write(f"{'='*80}\n")
                if df_sec.empty:
                    f.write("  No matches found.\n")
                    return
                hdr = (f"    {'Ligand_Name':<45} | "
                       f"{'Affinity':>10} | "
                       f"{'CNN_Pose':>10} | "
                       f"{'no_2og_binary':>13} | "
                       f"{'no_2og_pred_value':>18} | "
                       f"{'no_2og_ligand_iptm':>18}")
                f.write(hdr + "\n")
                f.write("    " + "-" * (len(hdr) - 4) + "\n")
                for _, row in df_sec.iterrows():
                    f.write(
                        f"    {row['Ligand_Name']:<45} | "
                        f"{row[affinity_col]:>10.2f} | "
                        f"{row[cnn_col]:>10.4f} | "
                        f"{row['no_2og_binary']:>13.6f} | "
                        f"{row['no_2og_pred_value']:>18.6f} | "
                        f"{row['no_2og_ligand_iptm']:>18.6f}\n"
                    )
                n_negative = int((df_sec['no_2og_pred_value'] < 0).sum())
                mean_binary   = df_sec['no_2og_binary'].mean()
                mean_pred_val = df_sec['no_2og_pred_value'].mean()
                f.write(f"\n  Total matched              : {len(df_sec)} compounds\n")
                f.write(f"  Negative no_2og_pred_value : {n_negative} / {len(df_sec)}"
                        f"  ({100 * n_negative / len(df_sec):.1f}%)\n")
                f.write(f"  Mean no_2og_binary         : {mean_binary:.6f}\n")
                f.write(f"  Mean no_2og_pred_value     : {mean_pred_val:.6f}\n")

            f.write(f"TOP-20% ({dataset_label}) – BOLTZ ENRICHMENT REPORT\n")
            write_section(enriched_affinity, "Top 20% by Vina Affinity – matched in boltz_dt")
            write_section(enriched_cnn,      "Top 20% by CNN Pose Score – matched in boltz_dt")

        print(f"Boltz-enriched table saved to: '{out_path}'")

    def run_boltz_enrichment_analysis(self):
        """Boltz enrichment for the NEW (BPP) dataset → top20_boltz_enriched.txt"""
        self._run_boltz_enrichment(
            txt_filename='top_20_new.txt',
            col_suffix='_new',
            dataset_label='NEW BPP DATASET',
            out_filename='top20_boltz_enriched.txt'
        )

    def run_boltz_enrichment_analysis_old(self):
        """Boltz enrichment for the OLD dataset → top20_boltz_enriched_old.txt"""
        self._run_boltz_enrichment(
            txt_filename='top_20_old.txt',
            col_suffix='_old',
            dataset_label='OLD DATASET',
            out_filename='top20_boltz_enriched_old.txt'
        )


    def run_cnn_affinity_enrichment(self):
        """New screening: Top 20% by CNN_Affinity from the raw BPP CSV,
        matched against boltz_dt.  Output: top20_cnn_affinity_boltz_enriched.txt"""
        if self.boltz_path is None:
            print("[WARN] boltz_path not set – skipping CNN_Affinity enrichment.")
            return
        if not os.path.exists(self.boltz_path):
            raise FileNotFoundError(f"boltz_dt.csv not found at '{self.boltz_path}'.")
        if self.df_bpp is None:
            raise RuntimeError("Call load_data() before run_cnn_affinity_enrichment().")

        # ------------------------------------------------------------------ #
        # 1. Pull Top 20% by CNN_Affinity (highest = best) from raw BPP data  #
        # ------------------------------------------------------------------ #
        df_bpp = self.df_bpp.copy()
        df_bpp['CNN_Affinity'] = pd.to_numeric(df_bpp['CNN_Affinity'], errors='coerce')
        df_bpp = df_bpp.dropna(subset=['CNN_Affinity'])

        top20_count = max(1, int(len(df_bpp) * 0.20))
        df_top = df_bpp.sort_values('CNN_Affinity', ascending=False).head(top20_count).copy()

        print(f"  [CNN_Affinity] {len(df_bpp)} valid rows in BPP dataset "
              f"→ top 20% = {top20_count} compounds")

        # ------------------------------------------------------------------ #
        # 2. Match against boltz_dt                                            #
        # ------------------------------------------------------------------ #
        boltz_cols = ['ligand_id', 'no_2og_binary', 'no_2og_pred_value', 'no_2og_ligand_iptm']
        df_boltz = pd.read_csv(self.boltz_path, usecols=boltz_cols)

        df_top['_lookup_id'] = df_top['Ligand_Name'].str.replace(
            r'_model_0_ligand$', '', regex=True
        )
        merged = df_top.merge(
            df_boltz,
            left_on='_lookup_id',
            right_on='ligand_id',
            how='inner'
        ).drop(columns=['_lookup_id', 'ligand_id'])

        # Preserve CNN_Affinity sort order after merge
        merged = merged.sort_values('CNN_Affinity', ascending=False)

        print(f"  [CNN_Affinity] {top20_count} top-20% entries "
              f"→ {len(merged)} matched in boltz_dt")

        # ------------------------------------------------------------------ #
        # 3. Write output report                                               #
        # ------------------------------------------------------------------ #
        out_filename = 'top20_cnn_affinity_boltz_enriched.txt'
        out_path = os.path.join(self.output_dir, out_filename)

        with open(out_path, 'w', encoding='utf-8') as f:
            f.write("TOP-20% by CNN_Affinity (NEW BPP DATASET) – BOLTZ ENRICHMENT REPORT\n")
            f.write(f"  Source file : {self.bpp_path}\n")
            f.write(f"  Total rows  : {len(df_bpp)}   |   Top 20% threshold: {top20_count} compounds\n")
            f.write(f"  Matched in boltz_dt: {len(merged)} compounds\n")
            f.write("\n")

            title = "Top 20% by CNN_Affinity – matched in boltz_dt"
            f.write(f"{'='*80}\n")
            f.write(f"  {title}\n")
            f.write(f"{'='*80}\n")

            if merged.empty:
                f.write("  No matches found.\n")
            else:
                # Determine which score columns are available
                vina_col = 'Vina_Affinity(kcal/mol)'
                cnn_pose_col = 'CNN_Pose_Score'
                cnn_aff_col  = 'CNN_Affinity'

                hdr = (f"    {'Ligand_Name':<45} | "
                       f"{'Vina_Affinity':>13} | "
                       f"{'CNN_Pose':>10} | "
                       f"{'CNN_Affinity':>12} | "
                       f"{'no_2og_binary':>13} | "
                       f"{'no_2og_pred_value':>18} | "
                       f"{'no_2og_ligand_iptm':>18}")
                f.write(hdr + "\n")
                f.write("    " + "-" * (len(hdr) - 4) + "\n")

                for _, row in merged.iterrows():
                    vina_str = f"{row[vina_col]:>13.2f}" if pd.notna(row.get(vina_col)) else f"{'N/A':>13}"
                    pose_str = f"{row[cnn_pose_col]:>10.4f}" if pd.notna(row.get(cnn_pose_col)) else f"{'N/A':>10}"
                    f.write(
                        f"    {row['Ligand_Name']:<45} | "
                        f"{vina_str} | "
                        f"{pose_str} | "
                        f"{row[cnn_aff_col]:>12.3f} | "
                        f"{row['no_2og_binary']:>13.6f} | "
                        f"{row['no_2og_pred_value']:>18.6f} | "
                        f"{row['no_2og_ligand_iptm']:>18.6f}\n"
                    )

                # Summary statistics
                n_negative = int((merged['no_2og_pred_value'] < 0).sum())
                mean_binary   = merged['no_2og_binary'].mean()
                mean_pred_val = merged['no_2og_pred_value'].mean()
                mean_cnn_aff  = merged[cnn_aff_col].mean()

                f.write(f"\n  Total matched              : {len(merged)} compounds\n")
                f.write(f"  Mean CNN_Affinity          : {mean_cnn_aff:.3f}\n")
                f.write(f"  Negative no_2og_pred_value : {n_negative} / {len(merged)}"
                        f"  ({100 * n_negative / len(merged):.1f}%)\n")
                f.write(f"  Mean no_2og_binary         : {mean_binary:.6f}\n")
                f.write(f"  Mean no_2og_pred_value     : {mean_pred_val:.6f}\n")

        print(f"CNN_Affinity-enriched table saved to: '{out_path}'")


    # Boltz no_2og_binary top-20% overlap analysis
    # ------------------------------------------------------------------
    def run_boltz_binary_top20_overlap(self):
        """Find top 20% of ALL boltz_dt rows by no_2og_binary (highest),
        then check which of those ligands appear in each of the three
        screening criteria independently:
          [Vina]  Top 20% by Vina Affinity      (from top20_boltz_enriched.txt, section 1)
          [CNN]   Top 20% by CNN Pose Score     (from top20_boltz_enriched.txt, section 2)
          [AFF]   Top 20% by CNN Affinity       (from top20_cnn_affinity_boltz_enriched.txt)
        Writes: boltz_binary_top20_overlap.txt
        """
        if self.boltz_path is None:
            print("[WARN] boltz_path not set – skipping binary overlap analysis.")
            return
        if not os.path.exists(self.boltz_path):
            raise FileNotFoundError(f"boltz_dt.csv not found at '{self.boltz_path}'.")

        # ── 1. Build top-20% of boltz_dt by no_2og_binary ────────────────────
        bcols = ['ligand_id', 'no_2og_binary', 'no_2og_pred_value', 'no_2og_ligand_iptm']
        df_boltz = pd.read_csv(self.boltz_path, usecols=bcols)
        df_boltz['no_2og_binary'] = pd.to_numeric(df_boltz['no_2og_binary'], errors='coerce')
        df_boltz = df_boltz.dropna(subset=['no_2og_binary'])

        top20_n = max(1, int(len(df_boltz) * 0.20))
        df_top  = df_boltz.sort_values('no_2og_binary', ascending=False).head(top20_n).copy()
        boltz_ids = set(df_top['ligand_id'].str.strip())

        print(f"  [binary overlap] {len(df_boltz)} boltz rows → top 20% = {top20_n} entries")

        # ── 2. Parse section-aware ligand ids from top20_boltz_enriched.txt ──
        def _parse_sections(filepath):
            """Returns dict: {'Vina Affinity': set(), 'CNN Pose Score': set()}
            by detecting the section headers inside the file."""
            sections = {'Vina Affinity': set(), 'CNN Pose Score': set()}
            current = None
            row_re  = re.compile(r'^\s+(P\S+)\s*\|')
            if not os.path.exists(filepath):
                print(f"  [WARN] '{filepath}' not found – skipping.")
                return sections
            with open(filepath, 'r', encoding='utf-8') as fh:
                for line in fh:
                    stripped = line.strip()
                    if 'Top 20% by Vina Affinity' in stripped:
                        current = 'Vina Affinity'
                        continue
                    if 'Top 20% by CNN Pose Score' in stripped:
                        current = 'CNN Pose Score'
                        continue
                    if current is None:
                        continue
                    m = row_re.match(line)
                    if m:
                        full = m.group(1).strip()
                        sections[current].add(re.sub(r'_model_0_ligand$', '', full))
            return sections

        def _parse_all_ids(filepath):
            """Returns flat set of all ligand ids found in a file."""
            ids = set()
            if not os.path.exists(filepath):
                print(f"  [WARN] '{filepath}' not found – skipping.")
                return ids
            row_re = re.compile(r'^\s+(P\S+)\s*\|')
            with open(filepath, 'r', encoding='utf-8') as fh:
                for line in fh:
                    m = row_re.match(line)
                    if m:
                        full = m.group(1).strip()
                        ids.add(re.sub(r'_model_0_ligand$', '', full))
            return ids

        enriched_path = os.path.join(self.output_dir, 'top20_boltz_enriched.txt')
        cnn_aff_path  = os.path.join(self.output_dir, 'top20_cnn_affinity_boltz_enriched.txt')

        sects       = _parse_sections(enriched_path)
        vina_ids    = sects['Vina Affinity']
        cnn_pose_ids = sects['CNN Pose Score']
        cnn_aff_ids  = _parse_all_ids(cnn_aff_path)

        print(f"  [binary overlap] Vina Affinity section    : {len(vina_ids)} ligands")
        print(f"  [binary overlap] CNN Pose Score section   : {len(cnn_pose_ids)} ligands")
        print(f"  [binary overlap] CNN Affinity file        : {len(cnn_aff_ids)} ligands")

        # ── 3. Compute overlaps for each criterion separately ─────────────────
        ov_vina     = boltz_ids & vina_ids       # boltz top-20% ∩ Vina Affinity top-20%
        ov_cnn_pose = boltz_ids & cnn_pose_ids   # boltz top-20% ∩ CNN Pose Score top-20%
        ov_cnn_aff  = boltz_ids & cnn_aff_ids    # boltz top-20% ∩ CNN Affinity top-20%

        # Combined views
        ov_all_three = ov_vina & ov_cnn_pose & ov_cnn_aff   # in all four lists
        ov_any       = ov_vina | ov_cnn_pose | ov_cnn_aff   # in at least one

        def _subset(id_set):
            return df_top[df_top['ligand_id'].isin(id_set)] \
                       .sort_values('no_2og_binary', ascending=False)

        # ── 4. Write report ───────────────────────────────────────────────────
        out_path = os.path.join(self.output_dir, 'boltz_binary_top20_overlap.txt')
        SEP  = '=' * 90
        SEP2 = '─' * 90

        def _write_table(fh, df, label, note=''):
            fh.write(f"\n{SEP}\n")
            fh.write(f"  {label}  ({len(df)} ligands){note}\n")
            fh.write(f"{SEP}\n")
            if df.empty:
                fh.write("  (none)\n")
                return
            hdr = (f"  {'ligand_id':<45} | "
                   f"{'no_2og_binary':>13} | "
                   f"{'no_2og_pred_value':>18} | "
                   f"{'no_2og_ligand_iptm':>18}")
            fh.write(hdr + "\n")
            fh.write("  " + "-" * (len(hdr) - 2) + "\n")
            for _, row in df.iterrows():
                fh.write(
                    f"  {row['ligand_id']:<45} | "
                    f"{row['no_2og_binary']:>13.6f} | "
                    f"{row['no_2og_pred_value']:>18.6f} | "
                    f"{row['no_2og_ligand_iptm']:>18.6f}\n"
                )

        with open(out_path, 'w', encoding='utf-8') as f:
            f.write("BOLTZ no_2og_binary TOP-20% — OVERLAP ANALYSIS\n")
            f.write(f"  boltz_dt total rows  : {len(df_boltz)}\n")
            f.write(f"  Top-20% threshold    : {top20_n} entries"
                    f" (no_2og_binary >= {df_top['no_2og_binary'].min():.6f})\n")
            f.write("\n  Three screening criteria checked independently:\n")
            f.write(f"    [Vina]  Top 20% by Vina Affinity    "
                    f"(top20_boltz_enriched.txt  §1) → {len(vina_ids)} ligands\n")
            f.write(f"    [CNN]   Top 20% by CNN Pose Score   "
                    f"(top20_boltz_enriched.txt  §2) → {len(cnn_pose_ids)} ligands\n")
            f.write(f"    [AFF]   Top 20% by CNN Affinity     "
                    f"(top20_cnn_affinity_boltz_enriched.txt) → {len(cnn_aff_ids)} ligands\n")

            # ── Summary table ─────────────────────────────────────────────
            f.write(f"\n{SEP2}\n")
            f.write(f"  OVERLAP SUMMARY\n")
            f.write(f"{SEP2}\n")
            f.write(f"  Boltz top-20% ∩ [Vina]              : {len(ov_vina):>4}\n")
            f.write(f"  Boltz top-20% ∩ [CNN Pose]          : {len(ov_cnn_pose):>4}\n")
            f.write(f"  Boltz top-20% ∩ [CNN Affinity]      : {len(ov_cnn_aff):>4}\n")
            f.write(f"  {SEP2}\n")
            f.write(f"  Boltz top-20% ∩ [Vina] ∩ [CNN Pose]: {len(ov_vina & ov_cnn_pose):>4}\n")
            f.write(f"  Boltz top-20% ∩ [Vina] ∩ [AFF]     : {len(ov_vina & ov_cnn_aff):>4}\n")
            f.write(f"  Boltz top-20% ∩ [CNN] ∩ [AFF]      : {len(ov_cnn_pose & ov_cnn_aff):>4}\n")
            f.write(f"  {SEP2}\n")
            f.write(f"  Boltz top-20% ∩ ALL THREE criteria  : {len(ov_all_three):>4}  ← strongest hits\n")
            f.write(f"  Boltz top-20% ∩ ANY criteria        : {len(ov_any):>4}\n")

            # ── Per-criterion tables ──────────────────────────────────────
            _write_table(f, _subset(ov_vina),
                         "Boltz top-20%  ∩  [Vina]  — Top 20% by Vina Affinity")
            _write_table(f, _subset(ov_cnn_pose),
                         "Boltz top-20%  ∩  [CNN Pose]  — Top 20% by CNN Pose Score")
            _write_table(f, _subset(ov_cnn_aff),
                         "Boltz top-20%  ∩  [AFF]  — Top 20% by CNN Affinity")
            _write_table(f, _subset(ov_all_three),
                         "Boltz top-20%  ∩  ALL THREE criteria  — highest confidence hits",
                         "  ← in every screening")

            # ── Plain name lists for quick lookup ─────────────────────────
            f.write(f"\n{SEP}\n")
            f.write(f"  PLAIN NAME LIST — Boltz top-20% ∩ [Vina]  ({len(ov_vina)} ligands)\n")
            f.write(f"{SEP}\n")
            for name in sorted(ov_vina):
                f.write(f"  {name}\n")

            f.write(f"\n{SEP}\n")
            f.write(f"  PLAIN NAME LIST — Boltz top-20% ∩ [CNN Pose]  ({len(ov_cnn_pose)} ligands)\n")
            f.write(f"{SEP}\n")
            for name in sorted(ov_cnn_pose):
                f.write(f"  {name}\n")

            f.write(f"\n{SEP}\n")
            f.write(f"  PLAIN NAME LIST — Boltz top-20% ∩ [CNN Affinity]  ({len(ov_cnn_aff)} ligands)\n")
            f.write(f"{SEP}\n")
            for name in sorted(ov_cnn_aff):
                f.write(f"  {name}\n")

            f.write(f"\n{SEP}\n")
            f.write(f"  PLAIN NAME LIST — in ALL THREE criteria  ({len(ov_all_three)} ligands)\n")
            f.write(f"{SEP}\n")
            if ov_all_three:
                for name in sorted(ov_all_three):
                    f.write(f"  {name}\n")
            else:
                f.write("  (none)\n")

            f.write(f"\n{SEP}\n")
            f.write(f"  PLAIN NAME LIST — in ANY criteria  ({len(ov_any)} ligands)\n")
            f.write(f"{SEP}\n")
            f.write("  Legend:  [V]=Vina  [C]=CNN Pose  [A]=CNN Affinity\n\n")
            for name in sorted(ov_any):
                tags = ""
                if name in ov_vina:     tags += " [V]"
                if name in ov_cnn_pose: tags += " [C]"
                if name in ov_cnn_aff:  tags += " [A]"
                f.write(f"  {name}{tags}\n")

        print(f"Boltz binary overlap report saved to: '{out_path}'")
        print(f"  → Boltz top-20% ∩ Vina Affinity    : {len(ov_vina)}")
        print(f"  → Boltz top-20% ∩ CNN Pose Score   : {len(ov_cnn_pose)}")
        print(f"  → Boltz top-20% ∩ CNN Affinity     : {len(ov_cnn_aff)}")
        print(f"  → In ALL three criteria             : {len(ov_all_three)}")
        print(f"  → In ANY criteria                   : {len(ov_any)}")


if __name__ == "__main__":
    # Define file paths dynamically
    base_dir  = r'c:\Users\Max\Desktop\Project\Comparing_docking_results'
    bpp_p     = os.path.join(base_dir, r'data\bpp_c\screen_summary_bpp.csv')
    old_p     = os.path.join(base_dir, r'data\old_c\screen_summary_old.csv')
    boltz_p   = os.path.join(base_dir, r'data\boltz_data\boltz_dt.csv')

    # Initialize and run the object-oriented comparer
    comparer = DockingComparer(
        bpp_path=bpp_p, old_path=old_p,
        output_dir=base_dir, boltz_path=boltz_p
    )
    comparer.load_data()
    comparer.generate_overall_summary()
    comparer.run_top_20_analysis()

    print("\n--- Boltz Enrichment Analysis (NEW) ---")
    comparer.run_boltz_enrichment_analysis()

    print("\n--- Boltz Enrichment Analysis (OLD) ---")
    comparer.run_boltz_enrichment_analysis_old()

    print("\n--- CNN_Affinity Top-20% Enrichment (NEW BPP) ---")
    comparer.run_cnn_affinity_enrichment()

    print("\n--- Boltz Binary Top-20% Overlap Analysis ---")
    comparer.run_boltz_binary_top20_overlap()
