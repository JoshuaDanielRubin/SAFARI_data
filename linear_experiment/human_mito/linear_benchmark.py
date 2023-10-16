import csv
import matplotlib.pyplot as plt

def calculate_metrics(TP, FP, TN, FN):
    sensitivity = TP / (TP + FN) if TP + FN != 0 else 0
    specificity = TN / (TN + FP) if TN + FP != 0 else 0
    precision = TP / (TP + FP) if TP + FP != 0 else 0
    accuracy = (TP + TN) / (TP + FP + TN + FN) if TP + FP + TN + FN != 0 else 0
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if precision + sensitivity != 0 else 0
    return sensitivity, specificity, precision, accuracy, f1_score

def report_metrics(metrics, label):
    print(f"{label}:")
    print(f"Sensitivity: {metrics[0]:.4f}")
    print(f"Specificity: {metrics[1]:.4f}")
    print(f"Precision: {metrics[2]:.4f}")
    print(f"Accuracy: {metrics[3]:.4f}")
    print(f"F1 Score: {metrics[4]:.4f}\n")

with open('alignment_stats.csv', 'r') as file:
    csv_reader = csv.DictReader(file)
    
    summary = {}
    
    for row in csv_reader:
        aligner = row["Aligner_Name"]
        if aligner not in summary:
            summary[aligner] = {
                "TP": 0, "FP": 0, "TN": 0, "FN": 0,
                "TP_MQ30": 0, "FP_MQ30": 0, "TN_MQ30": 0, "FN_MQ30": 0,
            }
        for key in summary[aligner]:
            summary[aligner][key] += int(row[key])
    
    aligner_order = ['safari', 'giraffe'] + [aligner for aligner in summary if aligner not in ['safari', 'giraffe']]
    metrics = ["Sensitivity", "Specificity", "Precision", "Accuracy", "F1 Score"]
    
    for metric_type, suffix in [("Overall", ""), ("MQ > 30", "_MQ30")]:
        print(f"------ {metric_type} ------\n")
        data_to_plot = {metric: [] for metric in metrics}
        
        for aligner in aligner_order:
            print(f"Aligner: {aligner}")
            results = calculate_metrics(summary[aligner]["TP" + suffix], summary[aligner]["FP" + suffix], summary[aligner]["TN" + suffix], summary[aligner]["FN" + suffix])
            report_metrics(results, metric_type)
            
            for metric in metrics:
                data_to_plot[metric].append(results[metrics.index(metric)])
        
        x = range(len(aligner_order))
        
        for metric in metrics:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.bar(x, data_to_plot[metric], label=metric, color=['b', 'r', 'g', 'y', 'c', 'm', 'k', 'orange', 'purple'])
            ax.set_xlabel('Aligners')
            ax.set_ylabel(metric)
            ax.set_title(f'{metric} ({metric_type}) by Aligner')
            ax.set_xticks(x)
            ax.set_xticklabels(aligner_order)
            ax.legend()
            plt.tight_layout()
            plt.savefig(f"linear_benchmark_{metric.replace(' ', '_').lower()}_{suffix}.png")
            plt.show()


