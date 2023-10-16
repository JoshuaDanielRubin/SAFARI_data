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
    
    # Sort the aligners to display Safari first and then Giraffe
    aligner_order = ['safari', 'giraffe'] + [aligner for aligner in summary if aligner not in ['safari', 'giraffe']]
    metrics = ["Sensitivity", "Specificity", "Precision", "Accuracy", "F1 Score"]
    
    for metric_type, suffix in [("Overall", ""), ("MQ > 30", "_MQ30")]:
        print(f"------ {metric_type} ------\n")
        data_to_plot = {metric: [] for metric in metrics}
        
        for aligner in aligner_order:
            print(f"Aligner: {aligner}")
            results = calculate_metrics(summary[aligner]["TP" + suffix], summary[aligner]["FP" + suffix], summary[aligner]["TN" + suffix], summary[aligner]["FN" + suffix])
            report_metrics(results, metric_type)
            
            for i, metric in enumerate(metrics):
                data_to_plot[metric].append(results[i])
        
        # Plotting
        x = range(len(aligner_order))
        width = 0.15
        fig, ax = plt.subplots(figsize=(12, 7))

        for i, metric in enumerate(metrics):
            ax.bar([pos + i * width for pos in x], data_to_plot[metric], width, label=metric)

        ax.set_xlabel('Aligners')
        ax.set_title(f'{metric_type} Metrics by Aligner')
        ax.set_xticks([pos + 2 * width for pos in x])
        ax.set_xticklabels(aligner_order)
        ax.legend()

        plt.tight_layout()
        plt.savefig("linear_benchmark.png")


