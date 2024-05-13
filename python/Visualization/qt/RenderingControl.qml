import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

ColumnLayout {
    id: columnLayout
    width: 100
    height: 200

    Switch {
        id: switch1
        text: qsTr("Switch")
        onCheckedChanged: {
            if(!checked) resize_details.start()
            else resize_details_open.start()
        }
    }

    Column {
        id: column

        Layout.fillWidth: true

        Rectangle {
            id: rectangle
            color: "#00fa0000"
            clip: true
            height: column.height
            width: column.width

            SequentialAnimation {
                id: resize_details
                NumberAnimation {
                    target: rectangle;
                    property: "height";
                    to: 0;
                    duration: 100
                }
            }
            SequentialAnimation {
                id: resize_details_open
                NumberAnimation {
                    target: rectangle;
                    property: "height";
                    to: 50;
                    duration: 100
                }
            }

            GridLayout {
                id: gridLayout
                anchors.fill: parent
                anchors.leftMargin: 5
                anchors.rightMargin: 5
                anchors.topMargin: 0
                anchors.bottomMargin: 0
                rows: 2
                columns: 3

                Label {
                    id: label1
                    text: qsTr("Label")
                    horizontalAlignment: Text.AlignHCenter
                    Layout.fillWidth: true
                }

                Slider {
                    id: slider
                    value: 0.5
                    Layout.fillWidth: true
                    onValueChanged: {
                        label3.text = value.toFixed(2)
                    }
                }


                Label {
                    id: label3
                    text: qsTr("Label")
                    horizontalAlignment: Text.AlignHCenter
                    Layout.fillWidth: true
                }

                Label {
                    id: label
                    text: qsTr("Label")
                    horizontalAlignment: Text.AlignHCenter
                    Layout.fillWidth: true
                }


                Slider {
                    id: slider1
                    value: 0.5
                    Layout.fillWidth: true
                    onValueChanged: {
                        label2.text = value.toFixed(2)
                    }
                }


                Label {
                    id: label2
                    text: qsTr("Label")
                    horizontalAlignment: Text.AlignHCenter
                    Layout.fillWidth: true
                }

                Item {
                    id: item2
                    Layout.fillHeight: true
                    Layout.fillWidth: true
                }
            }
        }
    }

    Switch {
        id: checkBox1
        opacity: 0.789
        text: qsTr("show triangle outline")
        display: AbstractButton.TextOnly
        Layout.fillWidth: true
        Layout.columnSpan: 1
        onCheckedChanged: function() {
            MainCtrl.toogleWireframe(checked)
        }
    }

    Switch {
        id: checkBox
        opacity: 0.789
        text: qsTr("white background")
        rotation: 0
        Layout.fillWidth: true
        Layout.columnSpan: 1
    }

    Item {
        id: item1
        Layout.fillHeight: true
        Layout.fillWidth: true
    }
}
