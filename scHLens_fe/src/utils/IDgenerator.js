import store from "@/store/index"

// //生产viewID
// export function generateViewID(){
//     let length = 6;
//     let characterList = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'];
//     let existViewIDs = store.state.dataList.map(item => item.ViewID)
//     let ID = '';
//     do{
//         ID = ''
//         for(let i = 0;i < length; i++){
//             ID += characterList[Math.floor(Math.random()*characterList.length)];
//         }
//     }
//     while(ID in existViewIDs)
//     return ID;
// }


