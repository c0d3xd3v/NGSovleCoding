import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Layouts 1.15

Item {
    MouseArea{
        anchors.fill: parent
        onClicked: {
            if(!root_item.focus){
                root_item.focus = true
            }else if(root_item.focus)
            {
            }
        }
        ColumnLayout{
            anchors.fill: parent
            /*
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
            }
            */
            GridLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                columns: 2
                Label{
                    text:"visualization"
                }
                ComboBox{
                    Layout.fillWidth: true
                    model: ["Source Mesh", "Result Mesh"]
                }
                Label{
                    text:"representation"
                }
                ComboBox{
                    Layout.fillWidth: true
                    model: ["Surface", "Wireframe", "Surface edges"]
                    onCurrentValueChanged: {
                        mesh_controller.setVisualization(currentValue);
                    }
                }
            }
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
            }
        }
    }
}
