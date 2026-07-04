import json
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

sns.set_theme(style="whitegrid")

try:
	with open("benchmark_results.json", "r") as f:
		data = json.load(f)
except Exception as e:
	print(f"Помилка читання файлу: {e}")
	exit(1)

methods_data = defaultdict(list)
for run in data:
	methods_data[run['method']].append(run)

for method_name, runs in methods_data.items():
	plt.figure(figsize=(10, 6))

	for run in runs:
		label = run['preconditioner']
		history = run['convergence_history']
		
		plt.semilogy(history, label=label, linewidth=1.8, marker='o', markersize=4, alpha=0.9)

	plt.title(f"Збіжність методу {method_name} з різними прекондиціонерами", fontsize=14, fontweight='bold', pad=15)
	plt.xlabel("Ітерації", fontsize=12)
	plt.ylabel("Норма нев'язки (log scale)", fontsize=12)
	plt.legend(loc="upper right", frameon=True, shadow=True)
	plt.tight_layout()

	filename_pdf = f"convergence_{method_name}.pdf"
	filename_png = f"convergence_{method_name}.png"
	plt.savefig(filename_pdf, format="pdf", bbox_inches='tight')
	plt.savefig(filename_png, dpi=300, bbox_inches='tight')
	print(f"📊 Графік для {method_name} збережено в '{filename_pdf}'")

	#plt.show()

print("\n" + "="*85)
print(f"{'Метод та прекондиціонер':<30} | {'Ітер.':<5} | {'Час солвера':<12} | {'Час прек.':<10} | {'Пам’ять (МБ)':<12}")
print("="*85)
for run in data:
	name = f"{run['method']} ({run['preconditioner']})"
	mem_total = run['memory_solve_mb'] + run['memory_precond_mb']
	mem_str = f"{mem_total:,.2f}" if mem_total > 0 else "-"
	print(f"{name:<30} | {run['iterations']:<5} | {run['time_solve']:<12.4f} | {run['time_precond']:<10.4f} | {mem_str:<12}")
print("="*85)
