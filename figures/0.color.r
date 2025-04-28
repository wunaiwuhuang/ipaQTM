colors_list <- list(
  "Palette 1"=c("#F56E55", "#DABA89", "#A4A4A4", "#4DDEC3", "#005873", "#AABA92", "#C27265", "#1977B1", "#8A8286", "#000070", "#E89F59", "#C6A6B3", "#000000", "#FFAE49", "#623D1A", "#C2CAC2", "#C95652", "#25494A", "#EF0079", "#73C2E9", "#7CBFAA", "#007A0B", "#76C3F7", "#00029E", "#60800A", "#8C9087", "#A4D0E3", "#D5B18E", "#873B00", "#53622A"),

  "Palette 2"=c("#0F4050", "#B06559", "#68855C", "#89243A", "#14696F", "#726F4C", "#ABACAC", "#CC4763", "#489003", "#D9AF6A", "#C8CEE4", "#D3D2D2", "#85C990", "#516B84", "#EB90BC", "#9D4097", "#D2E4A3", "#635477", "#3851A2", "#A16279", "#855D75", "#62666C", "#ED7823", "#EF4587"),

  "Palette 3"=c("#000000", "#9D9D9D", "#983E9A", "#CB9ECC", "#4776A6", "#ABC1D7", "#9E9C3D", "#CECD9E", "#973421", "#E5CCC7", "#F2589D", "#F8ABCE", "#58A44B", "#A2D8A2", "#4BA0F9", "#A5CFFC", "#F6A448", "#FAD1A3", "#489B9D", "#8B88B7", "#C0BFD8", "#B383B3", "#A3CDCE", "#D9D9D9"),

  "Palette 4"=c("#6EA0A0", "#74AB6D", "#EFA7AA", "#CC7276", "#2575B2", "#F08E40", "#F5BA6E", "#F07858", "#C3DDA5", "#DDA713", "#ED8760", "#349979", "#64BA9F", "#D62F81"),

  "Palette 5"=c("#006838", "#F56486", "#1F5577", "#00A651", "#EFCDCD", "#62A0BF", "#B5DF3C", "#908FB2", "#9A9999", "#A8DDB5", "#8B95F2", "#C6C6C5", "#BE1E2D", "#E2DCF0", "#EF4136", "#FAFAFA", "#E6E7E8")
  )
# 根据cancer名称生成调色板
generate_cancer_colors <- function(cancer) {

  # 将颜色和对应的标签展开为数据框
  color_df <- data.frame(
    Color = unlist(colors_list),
    Label = rep(names(colors_list), lengths(colors_list))
  )
  
  # 去重，只保留每种颜色第一次出现的标签
  unique_color_df <- color_df[!duplicated(color_df$Color), ]
  
  # 生成去重后的颜色向量和标签向量
  color_palette <- unique_color_df$Color
  color_labels <- unique_color_df$Label
  
  # 对每个cancer分配颜色
  cancer_colors <- setNames(color_palette[seq_along(cancer)], cancer)
  
  return(cancer_colors)
}