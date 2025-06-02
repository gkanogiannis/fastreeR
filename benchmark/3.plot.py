import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from PIL import Image

df = pd.read_csv("benchmark_logs.csv")
sns.set(style="whitegrid")
plt.rcParams.update({'figure.figsize': (12, 6), 'font.size': 12})
os.makedirs("plots", exist_ok=True)

df_fixed_samples = df[df['samples'] == 91]
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_fixed_samples, x='snps', y='runtime_seconds', hue='tool', style='threads', markers=True)
plt.title('Runtime vs SNPs (Samples = 91)')
plt.ylabel('Runtime (sec)')
plt.xlabel('Number of SNPs')
plot1_path = "plots/1.runtime_vs_snps_fixed_samples.png"
plt.tight_layout()
plt.savefig(plot1_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()

df_fixed_snps = df[df['snps'] == 2000000]
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_fixed_snps, x='samples', y='runtime_seconds', hue='tool', style='threads', markers=True)
plt.title('Runtime vs Samples (SNPs = 2M)')
plt.ylabel('Runtime (sec)')
plt.xlabel('Number of Samples')
plot2_path = "plots/2.runtime_vs_samples_fixed_snps.png"
plt.tight_layout()
plt.savefig(plot2_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()

df_fixed_config = df[(df['snps'] == 2000000) & (df['samples'] == 2000)]
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_fixed_config, x='threads', y='runtime_seconds', hue='tool', markers=True)
plt.title('Runtime vs Threads (SNPs = 2M, Samples = 2000)')
plt.ylabel('Runtime (sec)')
plt.xlabel('Threads')
plot3_path = "plots/3.runtime_vs_threads_fixed_config.png"
plt.tight_layout()
plt.savefig(plot3_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()

plt.figure(figsize=(10, 6))
sns.lineplot(data=df_fixed_samples, x='snps', y='memory_mb', hue='tool', style='threads', markers=True)
plt.title('Memory Usage vs SNPs (Samples = 91)')
plt.ylabel('Memory (MB)')
plt.xlabel('Number of SNPs')
plot4_path = "plots/4.memory_vs_snps_fixed_samples.png"
plt.tight_layout()
plt.savefig(plot4_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()

plt.figure(figsize=(10, 6))
sns.lineplot(data=df_fixed_snps, x='samples', y='memory_mb', hue='tool', style='threads', markers=True)
plt.title('Memory Usage vs Samples (SNPs = 2M)')
plt.ylabel('Memory (MB)')
plt.xlabel('Number of Samples')
plot5_path = "plots/5.memory_vs_samples_fixed_snps.png"
plt.tight_layout()
plt.savefig(plot5_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()

gradient_colors = [
    "#0570b0",   # deep blue
    "#3690c0",  # medium blue
    "#a6bddb",  # light blue
    "#fcbba1",  # light coral
    "#f7814a",  # orange-red 
    "#e41a1c"  # bright red
]
df_speedupA = df_fixed_snps.pivot_table(index=["snps", "samples", "threads"], columns="tool", values="runtime_seconds").reset_index()
df_speedupA["speedup"] = df_speedupA["vcf2dis"] / df_speedupA["fastreer"]
plt.figure(figsize=(10, 6))
palette = sns.color_palette("tab10") 
sns.lineplot(data=df_speedupA, x="samples", y="speedup", hue="threads", marker="o", palette=gradient_colors)
plt.title("Speedup of fastreeR (Java) over VCF2DIS (C++) (SNPs = 2M)")
plt.ylabel("Speedup (VCF2DIS / fastreeR)")
plt.xlabel("Number of Samples")
plt.axhline(1, color='gray', linestyle='--')
plot6A_path = "plots/6A.speedup.png"
plt.tight_layout()
plt.savefig(plot6A_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()
df_speedupB = df_fixed_samples.pivot_table(index=["snps", "samples", "threads"], columns="tool", values="runtime_seconds").reset_index()
df_speedupB["speedup"] = df_speedupB["vcf2dis"] / df_speedupB["fastreer"]
plt.figure(figsize=(10, 6))
palette = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
sns.lineplot(data=df_speedupB, x="snps", y="speedup", hue="threads", marker="o", palette=gradient_colors)
plt.title("Speedup of fastreeR (Java) over VCF2DIS (C++) (Samples = 91)")
plt.ylabel("Speedup (VCF2DIS / fastreeR)")
plt.xlabel("Number of SNPs")
plt.axhline(1, color='gray', linestyle='--')
plot6B_path = "plots/6B.speedup.png"
plt.tight_layout()
plt.savefig(plot6B_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()

subset_scale = df[df["snps"] == 2000000].copy()
subset_scale["config"] = subset_scale["samples"].astype(str) + " samples"
plt.figure(figsize=(12, 6))
sns.lineplot(data=subset_scale, x="threads", y="runtime_seconds", hue="config", style="tool", markers=True, dashes=False)
plt.xscale("log", base=2)
plt.title("Scalability vs Number of Threads (2M SNPs)")
plt.ylabel("Runtime (seconds)")
plot7_path = "plots/7.thread_efficieny.png"
plt.tight_layout()
plt.savefig(plot7_path, dpi=300, bbox_inches="tight")
#plt.show()
plt.close()


plot_files = [
    "plots/6B.speedup.png",
    "plots/6A.speedup.png",
    "plots/1.runtime_vs_snps_fixed_samples.png",
    "plots/2.runtime_vs_samples_fixed_snps.png"
]
images = [Image.open(p).convert("RGB") for p in plot_files]
w, h = images[0].size
images = [img.resize((w, h)) for img in images]
margin = 50
bg_color = (170, 170, 170)
canvas_width = w * 2 + margin * 3
canvas_height = h * 2 + margin * 3
combined_image = Image.new("RGB", (canvas_width, canvas_height), color=bg_color)
combined_image.paste(images[0], (margin, margin))
combined_image.paste(images[1], (w + margin*2, margin))
combined_image.paste(images[2], (margin, h + margin*2))
combined_image.paste(images[3], (w + margin*2, h + margin*2))
combined_image.save("plots/0.benchmark_combined.png", dpi=(300, 300))
