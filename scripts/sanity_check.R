
table(analysis1_df$response_4cat)

CR PD PR SD 
21 39 57 37 
> # Check recurrence counts by treatment group
+ analysis1_df %>%
+   filter(!is.na(recurrence_binary_sensitivity)) %>%
+   count(new_group, recurrence_binary_sensitivity) %>%
+   tidyr::pivot_wider(names_from = recurrence_binary_sensitivity, values_from = n, values_fill = 0)
# A tibble: 7 × 3
  new_group `No recurrence` Recurrence
  <chr>               <int>      <int>
1 CTx                    10         26
2 CTx + IO                9         13
3 CTx + RTx               2         13
4 IO                     11         15
5 Other                   5          1
6 Sx                     20         16
7 TKI                     4          9
> # Surgery pairs
+ length(surgery_pair_ids)
+ 
+ # CTx+IO pairs
+ paired_clinical %>%
+   filter(new_group == "CTx + IO") %>%
+   nrow()
[1] 13
[1] 4