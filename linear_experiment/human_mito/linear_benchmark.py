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
        damage_type = row["Damage_Type"]
        aligner = row["Aligner_Name"]
        
        if damage_type not in summary:
            summary[damage_type] = {}
        
        if aligner not in summary[damage_type]:
            summary[damage_type][aligner] = {
                "TP": 0, "FP": 0, "TN": 0, "FN": 0,
                "TP_MQ30": 0, "FP_MQ30": 0, "TN_MQ30": 0, "FN_MQ30": 0,
            }
        for key in summary[damage_type][aligner]:
            summary[damage_type][aligner][key] += int(row[key])

    aligner_order = ['safari', 'giraffe'] + [aligner for aligner in summary[next(iter(summary))].keys() if aligner not in ['safari', 'giraffe']]
    metrics = ["Sensitivity", "Specificity", "Precision", "Accuracy", "F1 Score"]
    
    for damage_type in summary:
        for metric_type, suffix in [("Overall", ""), ("MQ > 30", "_MQ30")]:
            data_to_plot = {metric: [] for metric in metrics}
            
            # Check which aligners have data for this damage type
            available_aligners = [aligner for aligner in aligner_order if aligner in summary[damage_type]]
            
            for aligner in available_aligners:
                results = calculate_metrics(summary[damage_type][aligner]["TP" + suffix], summary[damage_type][aligner]["FP" + suffix], summary[damage_type][aligner]["TN" + suffix], summary[damage_type][aligner]["FN" + suffix])
                
                # Print metrics to terminal
                print(f"Metrics for Damage Type: {damage_type}, Aligner: {aligner}, Metric Type: {metric_type}")
                for metric, value in zip(metrics, results):
                    print(f"{metric}: {value:.4f}")
                print("-" * 50)
                
                for metric in metrics:
                    data_to_plot[metric].append(results[metrics.index(metric)])
            
            x = np.arange(len(available_aligners))
            width = 0.15
            
            fig, ax = plt.subplots(figsize=(15, 7))
            
            for idx, metric in enumerate(metrics):
                ax.bar(x + idx*width, data_to_plot[metric], width, label=metric)

            ax.set_xlabel('Aligners')
            ax.set_ylabel('Score')
            ax.set_title(f'Metrics Comparison ({metric_type}) by Aligner for {damage_type}')
            ax.set_xticks(x + 2*width)
            ax.set_xticklabels(available_aligners)
            ax.legend()
            plt.ylim([0, 1])
            plt.tight_layout()

            if metric_type == "Overall":
                plt.savefig(f"{damage_type}_benchmark.png")
            else:
                plt.savefig(f"{damage_type}_benchmark_mq30.png")
            plt.show()

