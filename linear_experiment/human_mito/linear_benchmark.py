import csv
import numpy as np
import matplotlib.pyplot as plt

def calculate_metrics(TP, FP, TN, FN):
    sensitivity = TP / (TP + FN) if TP + FN != 0 else 0
    specificity = TN / (TN + FP) if TN + FP != 0 else 0
    precision = TP / (TP + FP) if TP + FP != 0 else 0
    accuracy = (TP + TN) / (TP + FP + TN + FN) if TP + FP + TN + FN != 0 else 0
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if precision + sensitivity != 0 else 0
    return sensitivity, specificity, precision, accuracy, f1_score

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
        data_to_plot = {metric: [] for metric in metrics}
        
        for aligner in aligner_order:
            results = calculate_metrics(summary[aligner]["TP" + suffix], summary[aligner]["FP" + suffix], summary[aligner]["TN" + suffix], summary[aligner]["FN" + suffix])
            
            for metric in metrics:
                data_to_plot[metric].append(results[metrics.index(metric)])
        
        x = np.arange(len(aligner_order))
        width = 0.15
        
        fig, ax = plt.subplots(figsize=(15, 7))
        
        for idx, metric in enumerate(metrics):
            ax.bar(x + idx*width, data_to_plot[metric], width, label=metric)

        ax.set_xlabel('Aligners')
        ax.set_ylabel('Score')
        ax.set_title(f'Metrics Comparison ({metric_type}) by Aligner')
        ax.set_xticks(x + 2*width)
        ax.set_xticklabels(aligner_order)
        ax.legend()
        plt.ylim([0, 1])
        plt.tight_layout()

        if metric_type == "Overall":
            plt.savefig("linear_benchmark.png")
        else:
            plt.savefig("linear_benchmark_mq30.png")
        plt.show()


